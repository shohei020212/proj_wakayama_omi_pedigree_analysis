#!/usr/bin/env nextflow

// DSL2
nextflow.enable.dsl = 2

// Module include statements
include { FASTQC as FASTQC_RAW    } from './modules/fastqc'
include { FASTQC as FASTQC_TRIM   } from './modules/fastqc'
include { TRIMMOMATIC             } from './modules/trimmomatic.nf'

workflow {

    // Check params
    if( !params.input ) error "Missing required parameter: --input (or set input in -params-file)"
    if( !params.fasta ) error "Missing required parameter: --fasta (or set fasta in -params-file)"

    def samplesheet = file(params.input)
    def ref_fasta   = file(params.fasta)

    // Create input channel from the contents of a CSV file
    channel
        .fromPath(samplesheet)
        .splitCsv(header:true)
        .map { row -> tuple(row.sample, file(row.fastq1), file(row.fastq2)) }
        .set { ch_reads }

    // Temporary diagnostics
    ch_reads.view()

    // Create a channel for the reference fasta
    channel.value(file(ref_fasta)).set { ch_ref }

    // Temporary diagnostics
    ch_ref.view()
    
    // 1) raw FastQC
    FASTQC_RAW(ch_reads)
    
    // 2) Trimming
    ch_adapters = params.trim_adapters ? channel.value(file(params.trim_adapters)) : channel.value(null)
    TRIMMOMATIC(ch_reads, ch_adapters)
    
    // 3) trimmed FastQC
    FASTQC_TRIM(TRIMMOMATIC.out.reads)
}
