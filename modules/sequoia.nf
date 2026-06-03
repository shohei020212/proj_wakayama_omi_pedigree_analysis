#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process SEQUOIA {
    // Save results to the specified output directory
    publishDir "${params.outdir}/sequoia", mode: 'copy'
    time = '24h'
    
    // Tag the process for easier monitoring
    tag "Running Sequoia on ${vcf_file.simpleName}"

    input:
    path vcf_file        // The VCF file
    path lh_file         // The LifeHistory file (or NO_FILE)

    output:
    path "Sequoia_Pedigree_${params.caller}.txt", emit: pedigree
    path "Sequoia_Diagnostic_Plots_${params.caller}.pdf"
    path "Sequoia_Summary_${params.caller}.txt"

    script:
    // Handle the optional file logic for the R script argument
    def lh_arg = lh_file.name != 'NO_FILE' ? "${lh_file}" : "NO_FILE"
    """
    sequoia_analysis.R ${vcf_file} ${lh_arg} ${params.caller}
    """
}