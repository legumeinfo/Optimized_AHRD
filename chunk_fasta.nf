#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_fasta = '/erdos/elavelle/ahrd_nextflow/glyma.Wm82.gnm6.ann1.PKSW.protein_primary.faa'
params.outdir = '/erdos/elavelle/ahrd_nextflow/Optimized_AHRD/'
params.chunksize = 500
params.database = ''

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
    python /erdos/elavelle/ahrd_nextflow/Optimized_AHRD/scripts/chunk_fasta.py ${input_fasta} ${out_dir} ${chunk_size}
    """
}

process alignChunks {
    tag "${chunk}"

    input: 
    path chunk
    val database

    output:
    path "diamond_out_${chunk}.txt"

    script:
    """
    export PATH=/home/elavelle/software:\$PATH

    //echo "Running: diamond blastp -d '${database}' -q '${chunk}' -o 'diamond_out_${chunk}.txt'"
    diamond blastp -d '${database}' -q '${chunk}' -o 'diamond_out_${chunk}.txt'
    """
}

workflow {
    //println "Database path: ${params.database}"
    input_fasta = file(params.input_fasta) 
    out_dir = params.outdir
    chunk_size = params.chunksize
    chunks = chunkFasta(input_fasta, out_dir, chunk_size)
    chunks 
        | flatten()
        | set { chunk_channel } 
    
    database_channel = Channel.value(params.database)
    alignChunks(chunk_channel, database_channel) 
}
