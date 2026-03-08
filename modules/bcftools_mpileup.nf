#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process BCFTOOLS_MPILEUP {
    tag "bcftools_joint_call"
    publishDir "${params.outdir}/bcftools", mode: 'copy'

    input:
    path bams
    path bais
    tuple path(ref_fasta), path(ref_index)
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
    -f ${ref_fasta} \
    ${bams} \
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
