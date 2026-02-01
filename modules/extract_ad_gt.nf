#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process EXTRACT_AD_GT {
    tag "extract_ad_gt"
    publishDir "${params.outdir}/genotype_qc", mode: 'copy'

    input:
    path filer_vcf

    output:
    path "ad_gt.tsv", emit: ad_gt_tsv

    script:
    """
    set -euo pipefail

    # Output columns: CHROM POS REF ALT SAMPLE GT AD
    bcftools query \
    -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%SAMPLE\\t%GT\\t%AD\\n]' \
    ${filer_vcf} > ad_gt.tsv
"""
}
