#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process PLOT_AD_SCATTER {
    tag "plot_ad_scatter"
    publishDir "${params.outdir}/genotype_qc", mode: 'copy'

    input:
    path ad_gt_tsv

    output:
    path "ad_gt_scatter.png"

    script:
    """
    # run bin/plot_ad_gt.py to generate scatter plot of AD
    plot_ad_gt.py --inputs ${ad_gt_tsv} --output ad_gt_scatter.png
    """
}
