#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process GATK_HC {
    tag "gatk_haplotype_caller"
    publishDir "${params.outdir}/gatk", mode: 'copy'
    time = '100h'

    input:
    path bams
    path bais
    tuple path(ref_fasta), path(ref_index)
    path regions_bed
    
    output:
    path "vars.vcf.gz", emit: vcf_gz

    script:
    def input_args = bams.collect { bam -> "-I ${bam}" }.join(" ")
    def regionOpt = regions_bed ? "-L ${regions_bed}" : ""
    """
    # Create sequence dictionary for the reference fasta
    gatk CreateSequenceDictionary -R ${ref_fasta} -O ${ref_fasta.baseName}.dict

    # Run GATK HaplotypeCaller
    gatk --java-options "-Xmx16g" HaplotypeCaller \
        --native-pair-hmm-threads ${task.cpus} \
        -R ${ref_fasta} \
        ${input_args} \
        -O vars.vcf.gz \
        ${regionOpt}
    """
}