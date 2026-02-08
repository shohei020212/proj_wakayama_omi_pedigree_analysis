#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process EXTRACT_AD_GT {
    tag "extract_ad_gt"
    publishDir "${params.outdir}/genotype_qc", mode: 'copy'

    input:
    path filer_vcf

    output:
    path "${params.caller}_ad_gt.tsv", emit: ad_gt_tsv

    script:
    """
    # output header
    echo -e "CHROM\tPOS\tREF\tALT\tSample\tGT\tAD" > ${params.caller}_ad_gt.tsv
    
    # append CHROM POS REF ALT SAMPLE GT AD for all samples
    bcftools query -f '[%CHROM\\t%POS\\t%REF\\t%ALT\\t%SAMPLE\\t%GT\\t%AD\\n]' ${filer_vcf} >> ${params.caller}_ad_gt.tsv
    """
}
