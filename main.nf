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
include { SAMTOOLS_VIEW_SORT      } from './modules/samtools'
include { MULTIQC_FLAGSTAT        } from './modules/multiqc_flagstat'

// Variant calling modules
include { BCFTOOLS_MPILEUP        } from './modules/bcftools_mpileup'
include { GATK_HC                 } from './modules/gatk_haplotypecaller'
include { FREEBAYES               } from './modules/freebayes'

include { VCFTOOLS_FILTER_QC      } from './modules/vcftools_filter_qc'
include { EXTRACT_AD_GT           } from './modules/extract_ad_gt'
include { PLOT_AD_SCATTER         } from './modules/plot_ad_scatter'
include { SEQUOIA                 } from './modules/sequoia'

workflow {

    // Check params
    if( !params.input ) error "Missing required parameter: --input (or set input in -params-file)"
    if( !params.fasta ) error "Missing required parameter: --fasta (or set fasta in -params-file)"
    if( !params.trim_adapters ) error "Missing required parameter: --adapter (or set adapter in -params-file)"
    
    def samplesheet = file(params.input)
    def ref_fasta   = file(params.fasta)
    def adapter     = file(params.trim_adapters)

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
    ch_adapters = channel.value(adapter)
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
    SAMTOOLS_VIEW_SORT(BWA.out.sam)
    
    // 6) Summarize flagstats with MultiQC
    ch_flagstats =  SAMTOOLS_VIEW_SORT.out.flagstats_filtered
        .mix(SAMTOOLS_VIEW_SORT.out.flagstats_unfiltered)
        .map { _sample, txt -> txt }
        .collect()
    
    MULTIQC_FLAGSTAT(ch_flagstats)

    // 7) Variant calling
    // Define channels for vcf output
    vcf_ch = channel.empty()
    // Collect BAM and BAI files
    ch_bam_list =  SAMTOOLS_VIEW_SORT.out.bam
        .map { _sample, bam -> bam }
        .collect()
    ch_bai_list =  SAMTOOLS_VIEW_SORT.out.bai
        .map { _sample, bai -> bai }
        .collect()

    ch_bam_list.view()
    ch_bai_list.view()
    
    // Variant calling step selection based on params.caller
    if (params.caller == 'bcftools') {
        
        // Optional regions BED file
        ch_regions = channel.value( file(params.regions, checkIfExists: true) )
        // Variant calling with BCFTOOLS mpileup + call
        vcf_ch = BCFTOOLS_MPILEUP(ch_bam_list, ch_bai_list, BWA_INDEX.out.ref_and_index, ch_regions).vcf_gz
        
    } else if (params.caller == 'gatk') {
        
        // Variant calling with GATK HaplotypeCaller
        vcf_ch = GATK_HC(ch_bam_list, ch_bai_list, BWA_INDEX.out.ref_and_index).vcf_gz
    
    } else if (params.caller == 'gatk') {
        
        // Variant calling with GATK HaplotypeCaller
        vcf_ch = GATK_HC(ch_bam_list, ch_bai_list, BWA_INDEX.out.ref_and_index).vcf_gz
    
    } else if (params.caller == 'freebayes') {
        
        // Variant calling with FREEBAYES
        vcf_ch = FREEBAYES(ch_bam_list, ch_bai_list, BWA_INDEX.out.ref_and_index).vcf_gz
    
    } else {

        error "Invalid variant caller specified: ${params.caller}. Please choose from 'bcftools', 'freebayes', or 'gatk'."
    
    }

    // 8) VCF filtering and QC metrics with VCFTOOLS
    //    Optional known_sites file
    ch_knownsites = channel.value( file(params.known_sites, checkIfExists: true) )
    VCFTOOLS_FILTER_QC(vcf_ch, ch_knownsites)

    // 9) Extract AD and GT from filtered VCF and plot AD scatter
    EXTRACT_AD_GT(VCFTOOLS_FILTER_QC.out.filtered_vcf)
    PLOT_AD_SCATTER(EXTRACT_AD_GT.out.ad_gt_tsv)

    // 10) Pedigree analysis with SEQUOIA
    ch_lhfile = channel.value( file(params.lh_file, checkIfExists: true) )
    SEQUOIA(VCFTOOLS_FILTER_QC.out.filtered_vcf, ch_lhfile)
}
