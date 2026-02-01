#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process BCFTOOLS_MPILEUP {
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
    path "vars.vcf.gz", emit: vcf_gz
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

    # Create renamed VCF to have simple sample names
    bcftools query -l vars.vcf.gz \
    | awk 'BEGIN{FS=OFS="\\t"}{
        old=\$0; new=\$0;
        sub(".*/","",new);              # delete path
        sub("\\\\.(bam|cram)\$","",new); # delete extension
        sub("([._-])sorted\$","",new);    # delete suffix "sorted"
        print old, new
        }' > sample_rename.tsv

    # reheader the VCF
    bcftools reheader -s sample_rename.tsv -o vars.renamed.vcf.gz vars.vcf.gz
    mv vars.renamed.vcf.gz vars.vcf.gz

    # index the renamed VCF file
    bcftools index -f vars.vcf.gz || true
    """
}
