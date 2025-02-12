#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_fasta = '/erdos/elavelle/ahrd_nextflow/glyma.Wm82.gnm6.ann1.PKSW.protein_primary.faa'
params.outdir = '/erdos/elavelle/ahrd_nextflow/'
params.chunksize = 500

process chunkFasta {
    tag { fasta }
    
    publishDir "${params.outdir}", mode: 'copy', overwrite: true 

    input:
    path input_fasta
    path out_dir
    val chunk_size

    output:
    stdout

    script:
    """
    python scripts/chunk_fasta.py ${input_fasta} ${out_dir} ${chunk_size}
    """
}

workflow {
    input_fasta = file( params.input_fasta ) 
    out_dir = params.outdir
    chunk_size = params.chunksize
    chunkFasta(input_fasta, out_dir, chunk_size)
}
