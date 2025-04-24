#!/usr/bin/env nextflow
nextflow.enable.dsl=2
import groovy.json.JsonOutput

params.chunksize = 500
params.maximum_processes = 128
params.threads = 4
params.databases = ''

simul_processes = params.maximum_processes / params.threads
diamond_processes = simul_processes / 2
interpro_processes = simul_processes / 2

assert params.input_fasta != ''
assert params.outdir != ''
assert params.interproscan_path != ''
assert params.diamond_path != ''
assert !params.databases.isEmpty(), "No database/s provided"

def check_params() {
    if( params.remove('help') ) {
        params.each{ k, v -> println "params.${k.padRight(25)} = ${v}" }
        exit 0
    }

    
    new File(params.outdir).mkdirs()
    if (!new File(params.outdir).isAbsolute()) {
        out_dir = new File(System.getenv('PWD'), params.outdir).getAbsolutePath()
        //def absoluteOutdir = new File(System.getenv('PWD'), params.outdir.toString()).getAbsolutePath()
        //params.outdir = absoluteOutdir
        println "INFO: outdir converted to absolute path: ${out_dir}"
    }
    new File(params.outdir).mkdirs()
    return out_dir
    // additional validation here
}

def load_database_csv(String csv_path) {
    def file = new File(csv_path)
    assert file.exists(), "CSV file '${csv_path}' does not exist."
    
    def databases_map = [:]
    file.eachLine { line, index ->
        if (line.trim() && !line.startsWith('#')) {
            def parts = line.split(',')
            assert parts.size() == 3 : "Line ${index + 1} in CSV must contain exactly 3 comma-separated values (dbName, dbPath (.fasta), dbIndex (.dmnd))"
            def (dbName, dbPath, dbIndex) = parts*.trim()
            databases_map[dbName] = [dbPath, dbIndex]
        }
    }
    return databases_map
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
    path myChunk

    output:
    path myChunk

    script:
    """
    sed -i 's/\\*/X/g' ${myChunk}
    """
}

process interproscan {
    tag "interproscan"
    
    input:
    tuple path(myChunk), path(out_dir), val(interproscan_path)

    output:
    path "${out_dir}/interproscan_out/${myChunk}.raw"

    maxForks interpro_processes

    script:
    """
    mkdir -p ${out_dir}/interproscan_out
    ${interproscan_path} -i $myChunk -f XML -o ${out_dir}/interproscan_out/${myChunk}.xml
    ${interproscan_path} -mode convert -f RAW -i ${out_dir}/interproscan_out/${myChunk}.xml -b ${out_dir}/interproscan_out/${myChunk}
    """
}

process concatenate_interproscan {
    tag "concatenate raw ips files"

    input:
    path raw_files
    path out_dir

    output:
    path "${out_dir}/interproscan_concatenated.raw"

    script:
    """
    cat ${raw_files} > ${out_dir}/interproscan_concatenated.raw
    """
}

process checkOutputExists {
    input:
    tuple val(db_name), val(db_path_str), val(db_index_str), path(outDir)

    output:
    tuple val(db_name), val(db_path_str), val(db_index_str), path("exists.txt")

    exec: //Need this for stdout message since NF hides stderr it in work dir log and also doesn't show echos or prints for processes that successfully finish -_-
    def checkFile = file("${outDir}/${db_name}_blasted.outfmt6.tsv").exists()
    if (checkFile) {
        log.warn("${outDir}/${db_name}_blasted.outfmt6.tsv already exists, will skip alignment to this database. YOU DO NOT WANT THIS IF YOU ARE USING A NEW QUERY FASTA IN AN OLD OUTDIR-either delete the outdir before rerunning, or provide another.")
    }
    script:
    """
    if [ -f "${outDir}/${db_name}_blasted.outfmt6.tsv" ]; then
        echo "true" >> exists.txt
    else
        echo "false" >> exists.txt
    fi
    """
}
//Currently not in use, functionality to check chunks (vs. concatenated file) to skip
process checkInterProExists {
    input:
    path out_dir
    path input_fasta

    output:
    tuple path(out_dir), path(input_fasta), path("IPSexists.txt")

    script:
    def chunk_name = input_fasta.name
    def out_file = "${out_dir}/interproscan_out/${chunk_name}.xml"
    """
    if [ -f "${out_file}" ]; then
        echo "true" > IPSexists.txt
    else
        echo "false" > IPSexists.txt
    fi
    """
}

process alignChunks {
    tag "diamond"

    input: 
    tuple path(myChunk), path(out_dir), val(database_name), val(database_index), val(diamond_path), val(database_path)

    output:
    tuple path("${out_dir}/diamond_out/${database_name}/${myChunk.baseName}.outfmt6.tsv"), val(database_name), val(database_path)
 
    maxForks diamond_processes    

    script:
    """
    mkdir -p ${out_dir}/diamond_out/${database_name}
    ${diamond_path} blastp -d ${database_index} -q ${myChunk} --outfmt 6 --threads ${params.threads} > ${out_dir}/diamond_out/${database_name}/${myChunk.getBaseName()}.outfmt6.tsv
    """
}

process concatenate_diamond_outputs {
    tag "concat diamonds"

    input:
    tuple path(diamond_in), path(out_dir), val(database_name), val(database_path)
    
    output:
    tuple path("${out_dir}/${database_name}_blasted.outfmt6.tsv"), val(database_name), val(database_path), emit: concatenated_files

    script:
    """
    cat ${diamond_in} > ${out_dir}/${database_name}_blasted.outfmt6.tsv
    """
}

process create_yaml {
    tag "Generate yaml"
    
    input:
    path fasta
    path out_dir

    output:
    path "${out_dir}/ahrd_config.yml"

    script:
    """
    python ${projectDir}/scripts/make_yaml.py ${fasta} ${out_dir} ${out_dir}/diamond_info.json
    
    cp ahrd_config.yml ${out_dir}/ahrd_config.yml
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
    
    ln -s /data/elavelle/databases/interpro.dtd
    java -jar /home/elavelle/software/AHRD/dist/ahrd.jar ahrd_config.yml
    cp ahrd_interpro_output.csv ${out_dir}/ahrd_output_file.csv
    """
}

workflow {
    println "${simul_processes}"
    input_fasta = file(params.input_fasta) 
    out_dir = check_params()

    chunk_size = params.chunksize
    def databases_map = load_database_csv(params.databases)
    log.info "LOADED DATABASES: ${databases_map}"

    chunks = chunkFasta(input_fasta, out_dir, chunk_size)
    clean_chunks = sanitize_fasta(chunks)
    clean_chunks 
        | flatten()
        | set { chunk_channel } 
    
    dbs_to_check = Channel.fromList(databases_map.collect { entry ->
        def dbName = entry.key
        def (dbPath, dbIndex) = entry.value
        tuple(dbName, file(dbPath), file(dbIndex), out_dir)
    })
    
    checkOutputExists(dbs_to_check)
        .branch { dbName, dbPath, dbIndex, exists_file ->
            exists: exists_file.text.trim() == "true"
            toProcess: exists_file.text.trim() == "false"
        }
        .set { existence_check }

    existence_check.toProcess
        .map { dbName, dbPath, dbIndex, exists_file ->
            tuple(dbName, dbPath, dbIndex)
        }
        .set { dbs_to_process }
    dbs_to_process.view { "DB TO PROCESS DEBUG: $it" }
 
    // Check if interproscan concatenated output already exists
    interproscan_result = file("${out_dir}/interproscan_concatenated.raw")
    if (interproscan_result.exists()) {
        log.warn("${out_dir}/interproscan_concatenated.raw already exists, will skip interproscan. YOU DO NOT WANT THIS IF YOU ARE USING A NEW QUERY FASTA IN AN OLD OUTDIR-either delete the outdir before rerunning, or provide another.")
        // Create a channel with the existing file
        interproscan_file_ch = Channel.value(file(interproscan_result))
    } else {
        // Run interproscan
        chunk_channel
            .map { chunk -> tuple(chunk, out_dir, params.interproscan_path) }
            .set { interproscan_input }
        //overwrite
        interproscan(interproscan_input)
            //.raw_file
            //.collect()
            .set { interproscan_output }

        interproscan_file_ch = concatenate_interproscan(interproscan_output, out_dir)
        //interproscan_file_ch = Channel.value(interproscan_result)
      }

    blast_input_channel = dbs_to_process
        .combine(chunk_channel)
        .map { dbName, dbPath, dbIndex, chunk ->
            tuple(chunk, out_dir, dbName, dbIndex, params.diamond_path, dbPath)
        }
    //blast_input_channel.view { "BLAST INPUT DEBUG: $it" }

    new_results = alignChunks(blast_input_channel)
        .map { result, dbName, dbPath ->
            tuple(result, dbName, dbPath)
        }

    existence_check.exists
        .map { dbName, dbPath, dbIndex, exists_file ->
            // Makes tuples with existing files that don't need alignment
            def existing_file = file("${out_dir}/${dbName}_blasted.outfmt6.tsv")
            tuple(existing_file, dbName, dbPath)
        }
        .set { existing_files }
    
    new_results // groups path and chunks by dbName
        .groupTuple(by: [1,2]) //dbName and dbPath
        .map { result, dbName, dbPath ->
            tuple(result, file(out_dir), dbName, dbPath)
        }
        .set { concat_input }

    concatenate_diamond_outputs(concat_input)
        .set { new_concatenated_files }
    
    //reunion
    //new_concatenated_files.view { println "FIRSTnu DEBUG: $it" }
    //existing_files.view { println "SECONDex DEBUG: $it" }
    def concatenated_files = existing_files.mix(new_concatenated_files)
   
    concatenated_files
        .map { blast_file, dbName, dbPath ->
            //def db_json = JsonOutput.toJson([dbName, blast_file.toString(), dbPath.toString()])
            [
                "dbName": dbName,
                "blast_file": blast_file.toString(),
                "dbPath": dbPath.toString()
            ]
        }
        .collect()
        .map { db_entry_list ->
            def db_json = groovy.json.JsonOutput.toJson(db_entry_list)
            def json_file = file("${out_dir}/diamond_info.json")       
            json_file.text = db_json
            tuple(
                input_fasta,
                json_file,
                out_dir
            )
        }
        // This is irrelevant, it returns null no matter what, in spite of the file being written
        .set { configInput }
     
    create_yaml(input_fasta, out_dir)
        .set { ahrd_config }

//    db_json_ch
//        .map { db_json ->
//            def dbs = new groovy.json.JsonSlurper().parseText(db_json)
//            def unsupported = dbs.find { !(it.dbName in ['uniprot_sprot', 'uniprot_trembl', 'uniref90']) }
//            if (unsupported) {
//                println "Unsupported database: ${unsupported.dbName} \nUser will need to manually add the appropriate regex lines in the AHRD .yaml"
//                System.exit(1)
//            }
//            db_json
//        }
//        .set { validated_db_json_ch }

    run_ahrd(out_dir, ahrd_config)
}
