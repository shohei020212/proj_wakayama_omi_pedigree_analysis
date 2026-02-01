#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process BCFTOOLS_MPILEUP_CALL {
    tag "bcftools_joint_call"
    publishDir "${params.outdir}/variants", mode: 'copy'

    cpus params.bcftools_threads ? params.bcftools_threads as int : 8
    memory params.bcftools_mem ? params.bcftools_mem : '16 GB'
    time params.bcftools_time ? params.bcftools_time : '8h'

    input:
    path bam_list
    path ref_fasta
    path regions_bed

    output:
    path "vars.vcf.gz"
    path "vars.vcf.gz.csi"
    
    script:
    def regionOpt = regions_bed ? "-R ${regions_bed}" : ""
    """
    set -euo pipefail

    # bcftools mpileup and call
    bcftools mpileup \
    --bam-list ${bam_list} \
    -f ${ref_fasta} \
    ${regionOpt} \
    -O b \
    -a AD,DP,SP \
    --threads ${task.cpus} \
    | bcftools call \
    -o vars.vcf.gz \
    -O z \
    -m \
    --threads ${task.cpus}

    # index the VCF file
    bcftools index -f vars.vcf.gz || true
    """
}
