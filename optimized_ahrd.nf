#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_fasta = ''
params.outdir = ''
params.chunksize = 500
params.database = ''
params.maximum_processes = 128
params.threads = 4 
simul_processes = params.maximum_processes / params.threads
diamond_processes = simul_processes / 2
interpro_processes = simul_processes / 2

check_params()

assert params.database != ''
assert params.input_fasta != ''
assert params.outdir != ''

def check_params() {
    if( params.remove('help') ) {
        params.each{ k, v -> println "params.${k.padRight(25)} = ${v}" }
        exit 0
    }
    // additional validation here
}

process chunkFasta {
    tag { fasta }
    
    publishDir "${params.outdir}", mode: 'copy', overwrite: true 

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
    
    output:
    path "interproscan_out/*.tsv"

    maxForks interpro_processes

    script:
    """
    mkdir -p "interproscan_out"

    /falafel/legumeinfo/sw/interproscan-5.72-103.0/interproscan.sh -i $myChunk -o interproscan_out/${myChunk}.tsv -f tsv
    """
}

process alignChunks {
    tag "${chunk}"

    input: 
    path chunk
    val database

    output:
    path "diamond_out_${chunk}.txt"
    
    maxForks diamond_processes    

    script:
    """
    export PATH=/home/elavelle/software:\$PATH

    diamond blastp -d '${database}' -q '${chunk}' -o 'diamond_out_${chunk}.txt' --threads '${params.threads}'
    """
}

workflow {
    //println "Database path: ${params.database}"
    println "${simul_processes}"
    input_fasta = file(params.input_fasta) 
    out_dir = params.outdir
    chunk_size = params.chunksize
    chunks = chunkFasta(input_fasta, out_dir, chunk_size)
    clean_chunks = sanitize_fasta(chunks)
    clean_chunks 
        | flatten()
        | set { chunk_channel } 
    
    database_channel = Channel.value(params.database)
  
    chunk_channel.set { chunk_files }
        interproscan(chunk_files)
        alignChunks(chunk_files, database_channel)
    
    //alignChunks(chunk_channel, database_channel) 
}
