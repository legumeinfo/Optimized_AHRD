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
    tag "${myChunk}"
    
    input:
    path myChunk
    path out_dir
    val interproscan_path 

    output:
    path "${out_dir}/interproscan_out/*.tsv"

    maxForks interpro_processes

    script:
    """
    mkdir -p ${out_dir}/interproscan_out
    ${interproscan_path} -i $myChunk -o ${out_dir}/interproscan_out/${myChunk}.tsv -f tsv
    """
}

process alignChunks {
    tag "${myChunk}"

    input: 
    path myChunk
    path out_dir
    val database
    val diamond_path

    output:
    path "${out_dir}/diamond_out/*.outfmt6.tsv"
    
    maxForks diamond_processes    

    shell:
    '''
    mkdir -p !{out_dir}/diamond_out
    !{diamond_path} blastp -d !{database} -q !{myChunk} --outfmt 6 --threads !{params.threads} > diamond_out_!{myChunk}.outfmt6.tsv
    '''
}

workflow {
    // !{out_dir}/diamond_out/diamond_out_${chunk_name}.outfmt6.tsv
    //println "Database path: ${params.database}"
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
  
    chunk_channel.set { chunk_files }
        interproscan(chunk_files, out_dir, params.interproscan_path)
        alignChunks(chunk_files, out_dir, database_channel, params.diamond_path)
    
    //alignChunks(chunk_channel, database_channel) 
}
