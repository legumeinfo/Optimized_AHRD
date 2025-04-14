#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_fasta = ''
params.outdir = ''
params.chunksize = 500
params.database = ''
params.maximum_processes = 128
params.threads = 4
params.interproscan_path = ''
params.diamond_path = '' 
simul_processes = params.maximum_processes / params.threads
diamond_processes = simul_processes / 2
interpro_processes = simul_processes / 2

// println file(params.diamond_path).getBaseName()

check_params()

assert params.database != ''
assert params.input_fasta != ''
assert params.outdir != ''
assert params.interproscan_path != ''
assert params.diamond_path != ''

def check_params() {
    if( params.remove('help') ) {
        params.each{ k, v -> println "params.${k.padRight(25)} = ${v}" }
        exit 0
    }
    if (!new File(params.outdir).isAbsolute()) {
        params.outdir = "$PWD/${params.outdir}"
        println "INFO: outdir converted to absolute path: ${params.outdir}"
    }
    // additional validation here
}

process chunkFasta {
    tag { fasta }

    input:
    path input_fasta
    path out_dir
    val chunk_size

    output:
    path "${out_dir}/chunks/*" 

    script:
    """
    python ${projectDir}/scripts/chunk_fasta.py ${input_fasta} ${out_dir} ${chunk_size}
    """
}

process sanitize_fasta {
    tag "fasta sanitization"

    input:
    path fasta_chunk

    output:
    path fasta_chunk

    script:
    """
    sed -i 's/\\*/X/g' ${fasta_chunk}
    """
}

process interproscan {
    tag "interproscan"
    
    input:
    path myChunk
    path out_dir
    val interproscan_path 

    output:
    path "${out_dir}/interproscan_out/${myChunk}.tsv"

    maxForks interpro_processes

    script:
    """
    mkdir -p ${out_dir}/interproscan_out
    ${interproscan_path} -i $myChunk -o ${out_dir}/interproscan_out/${myChunk}.tsv -f tsv
    """
}

process alignChunks {
    tag "diamond"

    input: 
    path myChunk
    path out_dir
    val database
    val diamond_path

    output:
    path "${out_dir}/diamond_out/diamond_out_${myChunk}.outfmt6.tsv", emit: chunk_out
    
    maxForks diamond_processes    

    shell:
    '''
    mkdir -p !{out_dir}/diamond_out
    !{diamond_path} blastp -d !{database} -q !{myChunk} --outfmt 6 --threads !{params.threads} > !{out_dir}/diamond_out/diamond_out_!{myChunk}.outfmt6.tsv
    '''
}

process concatenate_diamond_outputs {
    tag "concat diamonds"

    input:
    path diamond_cat_in
    path out_dir

    output:
    path "${out_dir}/diamond_concatenated.outfmt6.tsv"

    script:
    """
    cat ${diamond_cat_in} > ${out_dir}/diamond_concatenated.outfmt6.tsv
    """
}

process concatenate_interproscan_outputs {
    tag "concat interproscans"

    input:
    path interproscan_in
    path out_dir

    output:
    path "${out_dir}/interproscan_concatenated.tsv"

    script:
    """
    cat ${interproscan_in} > ${out_dir}/interproscan_concatenated.tsv
    """
}

process create_yaml {
    tag "Generate yaml"
    
    input:
    path proteins_fasta
    path diamond_file
    path interpro_result
    path database
    path out_dir

    output:
    path "${out_dir}/ahrd_config.yml"

    script:
    """
    python ${projectDir}/scripts/make_yaml.py ${proteins_fasta} ${diamond_file} ${interpro_result} ${database} ${out_dir}
    """
}

process run_ahrd{
    tag "Run ahrd"

    input:
    path out_dir
    path ahrd_config    

    output:
    path "ahrd_interpro_output.csv"

    script:
    """
    echo "" > blank1.txt
    echo "" > blank2.txt
    echo "" > blank3.txt

    java -jar /home/elavelle/software/AHRD/dist/ahrd.jar ahrd_config.yml
    cp ahrd_interpro_output.csv ${out_dir}/ahrd_output_file.csv
    """
}

workflow {
    println "${simul_processes}"
    input_fasta = file(params.input_fasta) 
    out_dir = file(params.outdir)
    out_dir.mkdirs()
    chunk_size = params.chunksize
    chunks = chunkFasta(input_fasta, out_dir, chunk_size)
    clean_chunks = sanitize_fasta(chunks)
    clean_chunks 
        | flatten()
        | set { chunk_channel } 
    
    database_channel = Channel.value(params.database)
  
    alignChunks(chunk_channel, out_dir, database_channel, params.diamond_path)
        .collect()
        .set { diamond_out_files }
    
    interproscan(chunk_channel, out_dir, params.interproscan_path)
        .collect()
        .set { interproscan_out_files }

    interproscan_out_files.view()

    diamond_cat = concatenate_diamond_outputs(diamond_out_files, out_dir)
    interproscan_cat = concatenate_interproscan_outputs(interproscan_out_files, out_dir)    
    
    create_yaml(input_fasta, diamond_cat, interproscan_cat, params.database, out_dir)
        .collect()
        .set { ahrd_config }

    run_ahrd(out_dir, ahrd_config)
}
