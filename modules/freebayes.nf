#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process FREEBAYES {
    tag "freebayes_joint_call"
    publishDir "${params.outdir}/freebayes", mode: 'copy'

    cpus params.freebayes_threads
    memory params.freebayes_mem
    time params.freebayes_time

    input:
    path bams
    path bais
    tuple path(ref_fasta), path(ref_index)

    output:
    path "vars.vcf.gz", emit: vcf_gz

    script:
    """
    freebayes -f ${ref_fasta} ${bams} > "vars.vcf"
    bgzip -@ ${task.cpus} vars.vcf
    """
}