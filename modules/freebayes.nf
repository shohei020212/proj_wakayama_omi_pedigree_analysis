#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process FREEBAYES {
    tag "freebayes_joint_call"
    publishDir "${params.outdir}/freebayes", mode: 'copy'
    time = '100h'
    memory = '20 GB'

    input:
    path bams
    path bais
    tuple path(ref_fasta), path(ref_index)
    path regions_bed

    output:
    path "vars.vcf.gz", emit: vcf_gz

    script:
    def regionOpt = regions_bed ? "--targets ${regions_bed}" : ""
    """
    freebayes -f ${ref_fasta} ${bams} ${regionOpt} > "vars.vcf"
    bgzip -@ ${task.cpus} vars.vcf
    """
}