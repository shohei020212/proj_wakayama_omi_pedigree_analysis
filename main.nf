#!/usr/bin/env nextflow

// DSL2
nextflow.enable.dsl = 2

// Module include statements
include { FASTQC as FASTQC_RAW    } from './modules/fastqc'
include { FASTQC as FASTQC_TRIM   } from './modules/fastqc'
include { TRIMMOMATIC             } from './modules/trimmomatic'
include { MULTIQC                 } from './modules/multiqc'
include { BWA_INDEX               } from './modules/bwa_index'
include { BWA                     } from './modules/bwa'
include { SAMTOOLS                } from './modules/samtools'
include { MULTIQC_FLAGSTAT        } from './modules/multiqc_flagstat'

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
    
    // Create a channel for the reference fasta
    channel.value(file(ref_fasta)).set { ch_ref }

    // Temporary diagnostics
    ch_reads.view()
    ch_ref.view()

    // 1) raw FastQC before trimming
    FASTQC_RAW(ch_reads)
    
    // 2) Trimming with Trimmomatic
    ch_adapters = params.trim_adapters ? channel.value(file(params.trim_adapters)) : channel.value(null)
    TRIMMOMATIC(ch_reads, ch_adapters)
    
    // 3) trimmed FastQC with Trimmomatic output
    FASTQC_TRIM(TRIMMOMATIC.out.reads)
    
    // 4) Summarize raw QC and trimmed QC with MultiQC
    //    Mix FASTQC output fies(zip/html）
    ch_fastqc_files = FASTQC_RAW.out
        .mix(FASTQC_TRIM.out)
        .map { _sample, zip, html -> [ zip, html ] }
        .flatten()
        
    MULTIQC(ch_fastqc_files.collect())

    // 4) Alignment with BWA-MEM
    BWA_INDEX(ch_ref)
    BWA(TRIMMOMATIC.out.reads, BWA_INDEX.out.ref_and_index)

    // 5) SAM to sorted BAM with SAMTOOLS
    SAMTOOLS(BWA.out.sam)
    
    // 6) Summarize flagstats with MultiQC
    ch_flagstats =  SAMTOOLS.out.flagstats_filtered
        .mix(SAMTOOLS.out.flagstats_unfiltered)
        .map { _sample, txt -> txt }
        .collect()
    
    MULTIQC_FLAGSTAT(ch_flagstats)
    
    // 

    // End of workflow

}
