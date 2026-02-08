#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process GATK_HC {
    tag "gatk_haplotype_caller"
    publishDir "${params.outdir}/gatk", mode: 'copy'

    cpus params.gatk_threads
    memory params.gatk_mem
    time params.gatk_time

    input:
    path bams
    path bais
    tuple path(ref_fasta), path(ref_index)

    output:
    path "vars.vcf.gz", emit: vcf_gz

    script:
    def input_args = bams.collect { bam -> "-I ${bam}" }.join(" ")
    """
    # Create sequence dictionary for the reference fasta
    gatk CreateSequenceDictionary -R ${ref_fasta} -O ${ref_fasta.baseName}.dict

    # Run GATK HaplotypeCaller
    gatk HaplotypeCaller \
        --native-pair-hmm-threads ${task.cpus} \
        -R ${ref_fasta} \
        ${input_args} \
        -O vars.vcf.gz
    """
}