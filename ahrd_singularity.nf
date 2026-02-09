#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import groovy.json.JsonOutput

params.chunksize = 500
params.maximum_processes = 128
params.threads = 4
params.databases = ''
params.gaf = ''
if( !params.containsKey('desc_blacklist') ) {
    params.desc_blacklist = '/resources/blacklist_descline.txt'
}
if( !params.containsKey('token_blacklist') ) {
    params.token_blacklist = '/resources/blacklist_token.txt'
}

simul_processes = params.maximum_processes / params.threads
diamond_processes = simul_processes / 2
interpro_processes = simul_processes / 2

assert params.input_fasta != ''
assert params.outdir != ''
assert !params.databases.isEmpty(), "No database/s provided"
assert params.gaf != ''

def check_params() {
    if( params.remove('help') ) {
        params.each{ k, v -> println "params.${k.padRight(25)} = ${v}" }
        exit 0
    }

    if (!new File(params.outdir).isAbsolute()) {
        out_dir = new File(System.getenv('PWD'), params.outdir).getAbsolutePath()
        //def absoluteOutdir = new File(System.getenv('PWD'), params.outdir.toString()).getAbsolutePath()
        //params.outdir = absoluteOutdir
        println "INFO: outdir converted to absolute path: ${out_dir}"
    }
    new File(out_dir).mkdirs()
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
    path "${out_dir}/chunks/*fasta" 

    script:
    """
    python3 /bin/scripts/chunk_fasta.py ${input_fasta} ${out_dir} ${chunk_size} 
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
    sed -i '/^>/! {s/\\*/X/g; s/\\./X/g}' ${myChunk}
    """
}

process interproscan {
    tag "interproscan"
    
    input:
    tuple path(myChunk), path(out_dir)

    output:
    path "${out_dir}/interproscan_out/${myChunk}.raw"

    maxForks interpro_processes

    script:
    """
    set -x
    #source /environment

    mkdir -p ${out_dir}/interproscan_out
    mkdir -p ${out_dir}/interproscan_tmp
    /opt/interproscan-5.72-103.0/interproscan.sh --disable-precalc -i $myChunk -f XML --goterms -o ${out_dir}/interproscan_out/${myChunk}.xml
    if [ -f "${out_dir}/interproscan_out/${myChunk}.raw" ]; then
        rm "${out_dir}/interproscan_out/${myChunk}.raw"
    fi
    /opt/interproscan-5.72-103.0/interproscan.sh -mode convert -f RAW -i ${out_dir}/interproscan_out/${myChunk}.xml -b ${out_dir}/interproscan_out/${myChunk}

    if [ ! -f "${out_dir}/interproscan_out/${myChunk}.raw" ]; then
        echo "ERROR: No hits to create IPS output. Add more .fasta sequences databases. Deleting interproscan_out directory."
        rm -r ${out_dir}/interproscan_out
        exit 1
    fi 
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
    tuple path(myChunk), path(out_dir), val(database_name), val(database_index), val(database_path)

    output:
    tuple path("${out_dir}/diamond_out/${database_name}/${myChunk.baseName}.outfmt6.tsv"), val(database_name), val(database_path)
 
    maxForks diamond_processes

    script:
    """
    mkdir -p ${out_dir}/diamond_out/${database_name}
    diamond blastp -d ${database_index} -q ${myChunk} --outfmt 6 --threads ${params.threads} > ${out_dir}/diamond_out/${database_name}/${myChunk.getBaseName()}.outfmt6.tsv
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
    tag "Generate yaml per fasta chunk"
    
    input:
    tuple path(fasta), path(json_chunk)
    path out_dir
    path gaf
    path ips_xml
    path desc_blacklist
    path token_blacklist
    path IPSdummy 

    output:
    path "${out_dir}/chunks/${json_chunk.baseName}_ahrd_config.yml"
    
    script:
    """
    python3 /bin/scripts/make_yaml.py ${fasta} ${out_dir} ${json_chunk} ${gaf} ${ips_xml} ${desc_blacklist} ${token_blacklist}

    mv ahrd_config.yml ${out_dir}/chunks/${json_chunk.baseName}_ahrd_config.yml
    """
}

process run_ahrd{
    tag "Run ahrd"

    input:
    each path(out_dir)
    path ahrd_config

    output:
    path "${out_dir}/chunks/${ahrd_config.baseName}_ahrd_output.csv", emit: ahrd_out

    script:
    """
    echo "" > blank1.txt
    ln -s /resources/interpro.xml
    ln -s /resources/interpro.dtd
    java -jar /bin/scripts/ahrd.jar ${out_dir}/chunks/${ahrd_config}
    cp ahrd_interpro_output.csv ${out_dir}/chunks/${ahrd_config.baseName}_ahrd_output.csv
    """
}

process merge_ahrd {
    tag "merge ahrd configs"

    input:
    path out_dir
    path configs

    output:
    path "${out_dir}/ahrd_merged.csv", emit: ahrd_out

    script:
    """
    # Grabs header
    grep -m 1 '^Protein-Accession' ${configs[0]} > merged.csv
    # Ignores empty lines and duplicate headers to append
    tail -n +2 -q ${configs.join(' ')} | grep -v '^Protein-Accession' | grep -v '^\$' >> merged.csv   
    #head -n 1 ${configs[0]} > merged.csv
    #tail -n +2 -q ${configs.join(' ')} >> merged.csv
    mv merged.csv ${out_dir}/ahrd_merged.csv
    """
}

process appendGOids {
    tag "tack on GO ids associated with queries via IPS"
    
    input:
    path ahrd_output
    path interpro_raw
    path out_dir

    output:
    path "${out_dir}/ahrd_with_go_ids.tsv"

    script:
    """
    #Debug - list all inputs
    echo "AHRD output file: ${ahrd_output}"
    echo "InterPro file: ${interpro_raw}" 
    echo "Output dir: ${out_dir}"
    
    perl /bin/scripts/ipr2go.pl $interpro_raw > go_map.tsv
    cp go_map.tsv ${out_dir}/
    echo "Perl script completed, go_map.tsv size:"
    ls -la go_map.tsv
    head -5 go_map.tsv

    awk 'BEGIN {FS=OFS="\t"}
        FNR==NR {
            key = \$1
            val = \$2 " (" \$3 ")"
            if (key in map) {
                map[key] = map[key] ", " val;
            } else {
                map[key] = val;
            }
            next;
    }
    FNR==1 {
        print \$0, "Query_GO_terms";  # Add header for GO terms
        next;
    }
    FNR > 3 {
        go_terms = (\$1 in map) ? map[\$1] : "";
        print \$0, go_terms;
    }' go_map.tsv "$ahrd_output" > "${out_dir}/ahrd_with_go_ids.tsv"
    echo "Final output size:"

    ls -la ${out_dir}/ahrd_with_go_ids.tsv

    # Fail if output is empty
    if [ ! -s "${out_dir}/ahrd_with_go_ids.tsv" ]; then
        echo "ERROR: Output file is empty! Check to see if the interproscan_concatenated.raw file has only NULL values in the last column."
        exit 1
        fi   
    """
}

go_obo_file="/resources/go.obo"

process create_GO_lookup {
    tag "Create lookup table from hardcoded .obo"

    input:
    //path go_obo_file
    path out_dir

    output:
    path "go_term_lookup.tsv"

    script:
    """
    echo "GO lookup entered using local file"
    # Download the .obo
    # curl -s --connect-timeout 60 --max-time 600 https://current.geneontology.org/ontology/go.obo -o /tmp/go.obo 2>curl_error.log || { echo "ERROR: Failed to download go.obo"
    # parse to make lookup table 
    awk '
    BEGIN { 
        FS=" ";
        print "GO_ID\\tDescription"; 
    }
    /^\\[Term\\]/ { term=1; id=""; name=""; next }
    term==1 && /^id: GO:/ { id=\$2; next }
    term==1 && /^name:/ { 
        name=substr(\$0, 7); 
        if (id != "" && name != "") {
            print id "\\t" name;
            id=""; name=""; term=0;
        }
    }
    ' ${go_obo_file} > go_term_lookup.tsv
    ls -ltrh go_term_lookup.tsv
    cp go_term_lookup.tsv ${out_dir}/go_term_lookup.tsv
    """
}

process description_for_hits {
    tag "Insert descriptions for hit GO terms"

    input:
    path ahrd_output_with_go
    path go_lookup
    path out_dir

    output:
    path "${out_dir}/ahrd_with_all_descriptions.tsv"

    script:
    """
    echo "input files: "
    ls -la ${ahrd_output_with_go} ${go_lookup}

    echo "input file true sizes: "
    ls -la `readlink ${ahrd_output_with_go}`
    ls -la `readlink ${go_lookup}`

    python3 /bin/scripts/add_hit_descriptors.py ${ahrd_output_with_go} ${go_lookup} ${out_dir}   
    """

    }

process write_json {
    publishDir "${params.outdir}/chunks", mode: 'copy'

    input:
    tuple val(chunkId), val(db_entries)
    
    output:
    path "diamond_info_${chunkId}.json"
    
    script:
    """
    cat > diamond_info_${chunkId}.json << 'EOF'
    ${groovy.json.JsonOutput.toJson(db_entries)}
    EOF
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
        interproscan_output = Channel.value(file(interproscan_result))
    } else {
        // Run interproscan
        chunk_channel
            .map { chunk -> tuple(chunk, out_dir) }
            .set { interproscan_input }
        //overwrite
        interproscan(interproscan_input)
            .collect()
            .set { interproscan_output }

        interproscan_file_ch = concatenate_interproscan(interproscan_output, out_dir)
      }

    blast_input_channel = dbs_to_process
        .combine(chunk_channel)
        .map { dbName, dbPath, dbIndex, chunk ->
            tuple(chunk, out_dir, dbName, dbIndex, dbPath)
        }

    new_files = alignChunks(blast_input_channel)
        .map { result, dbName, dbPath ->
            tuple(result, dbName, dbPath)
        }
    
    existence_check.exists
        .map { dbName, dbPath, dbIndex, exists_file ->
            // Makes tuples with existing files that don't need alignment
            tuple(dbName, dbPath)
        }
        .set { existing_files }
    
    new_files // groups path and chunks by dbName
        .groupTuple(by: [1,2]) //dbName and dbPath
        .map { result, dbName, dbPath ->
            tuple(result, file(out_dir), dbName, dbPath)
        }
        .set { concat_input }
    
    // Not strictly part of this pipeline, but creates files for the check if out_dir is reused
    concatenate_diamond_outputs(concat_input)
        //set {new_concatenated_files}
    //reunion
    //def concatenated_files = existing_files.mix(new_concatenated_files)
 
    //concatenated_files
        //.map { blast, db, path -> tuple(db, path) }
        //.collect()
        //.set { db_paths_list }

    // THIS WILL ONLY RUN ON EXISTING DATABASES (IT DOES NOT WAIT OR RUN IN SEQUENCE)
    Channel.fromPath("${out_dir}/diamond_out/*/*.outfmt6.tsv")
        .map { chunkBlast ->
            def dbName = chunkBlast.parent.name
            def chunkId = chunkBlast.name.split('\\.')[0]
            tuple(dbName, chunkId, file(chunkBlast))
        }
        .set { existsPart2 }
    
    existing_files
        .join(existsPart2)
        .map { dbNames, dbPath, chunkId, chunkBlast ->
            tuple(chunkId, chunkBlast, dbNames, dbPath)
        //println "DEBUG - existDBpaths: ${dbPath}"
        }
        .set { existsFull }

    new_files.map { chunkBlast, dbNames, dbPath ->
            chunkBlast = file(chunkBlast)
            def chunkId = chunkBlast.name.split('\\.')[0]
            tuple(chunkId, chunkBlast, dbNames, dbPath)
        }
        .set {just_new}

    existsFull
        .mix(just_new)
        .map { chunkId, chunkBlast, dbNames, dbPaths ->
         tuple(
            chunkId,
            chunkBlast,
            dbNames,
            dbPaths
        )
        }
        .set { all }
 
        all 
        // "compresses" list: e.g.:(chunkId, [blastFiles], [dbNames], dbName1, dbPath1, etc.)
        .groupTuple(by: 0)

        .map { args ->
            //println "DEBUG - full args: ${args}"
            def chunkId = args[0]
            def blastFiles = args[1]
            def dbNames = args[2]
            def dbPath = args[3]
   
            def db_entries = (0..<dbNames.size()).collect { i ->
 	         def dbName = dbNames[i]
	         ["dbName": dbName, "blast_file": "${out_dir}/diamond_out/${dbName}/${blastFiles[i].name}", "dbPath": dbPath[i]]
                }
                def json_file = file("${out_dir}/chunks/diamond_info_${chunkId}.json")
                def json_content = "[" + db_entries.collect { entry ->
                // MANUAL creation of .json, because the groovy function causes stack overflow for some reason
                """{"dbName":"${entry.dbName}","blast_file":"${entry['blast_file']}","dbPath":"${entry.dbPath}"}"""
                }.join(",") + "]"
                
                json_file.text = json_content
                def chunk_fasta = file("${out_dir}/chunks/${chunkId}.fasta")
                tuple(chunk_fasta, json_file)
        }
        .set {yaml_input_channel}

    create_yaml(yaml_input_channel, out_dir, file(params.gaf).toAbsolutePath().toString(), file("/resources/interpro.xml"), params.desc_blacklist, params.token_blacklist, interproscan_output)
        .set { ahrd_configs }

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
    
    // solely because nf refuses to use out_dir in a definition again, duplicate it
    def out_dir2 = check_params()
    def ahrd_output = file("${out_dir2}/").list().any { it.endsWith('hrd_merged.csv') } 
    //def ahrd_output = file("${out_dir2}/ahrd_merged.csv").size() > 0

    if (!ahrd_output) {
        run_ahrd(out_dir, ahrd_configs)
        .collect()
        .set { ahrd_chunks }

        merge_ahrd(out_dir, ahrd_chunks)
        .set { ahrd_output_ch }
    }
    else{
        log.warn("Merged AHRD result file already exists in out dir. Delete if you wish to regenerate it.")
        Channel.fromPath("${out_dir}/ahrd_merged.csv")
            .set { ahrd_output_ch }
    }

    appendGOids(ahrd_output_ch, interproscan_output, out_dir)
        .set { appendedOut }

    create_GO_lookup(out_dir)
        .set { go_lookup }

    description_for_hits(appendedOut, go_lookup, out_dir)
}
