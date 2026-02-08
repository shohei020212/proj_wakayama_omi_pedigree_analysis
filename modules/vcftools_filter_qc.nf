#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process VCFTOOLS_FILTER_QC {
    tag "vcftools"
    publishDir "${params.outdir}/vcftools", mode: 'copy'

    input:
    path vcf_gz
    path known_sites

    output:
    // filtered VCF
    path "snps.filtered.recode.vcf"
    path "snps.filtered.rm_miss_indv.recode.vcf"
    path "snps.filtered.rm_miss_indv.site.recode.vcf", emit: filtered_vcf

    // QC outputs
    path "snps.filtered.recode.frq"
    path "snps.filtered.recode.idepth"
    path "snps.filtered.recode.ldepth.mean"
    path "snps.filtered.recode.lqual"
    path "snps.filtered.recode.imiss"
    path "snps.filtered.recode.lmiss"
    path "remove_indv.txt" 

    script:
    def posOpt = (known_sites.name != 'NO_FILE') ? "--positions ${known_sites}" : ""
    """
    set -euo pipefail
    
    # -------------------------------------------------------
    # 1) SNP filtering
    # -------------------------------------------------------
    vcftools \
    --gzvcf ${vcf_gz} \
    ${posOpt} \
    --remove-indels \
    --min-alleles ${params.vcftools_min_alleles} \
    --max-alleles ${params.vcftools_max_alleles} \
    --min-meanDP  ${params.vcftools_min_meanDP} \
    --recode \
    --out snps.filtered
    # -------------------------------------------------------
    # 2) QC metrics
    # -------------------------------------------------------
    vcftools \
    --vcf snps.filtered.recode.vcf \
    --freq2 \
    --out snps.filtered.recode

    vcftools \
    --vcf snps.filtered.recode.vcf \
    --depth \
    --out snps.filtered.recode

    vcftools \
    --vcf snps.filtered.recode.vcf \
    --site-mean-depth \
    --out snps.filtered.recode

    vcftools \
    --vcf snps.filtered.recode.vcf \
    --site-quality \
    --out snps.filtered.recode

    vcftools \
    --vcf snps.filtered.recode.vcf \
    --missing-indv \
    --out snps.filtered.recode

    vcftools \
    --vcf snps.filtered.recode.vcf \
    --missing-site \
    --out snps.filtered.recode
    # -------------------------------------------------------
    # 3) Create a list of removed individuals due to missingness
    # -------------------------------------------------------
    awk '\$5 > ${params.vcftools_ind_max_miss} {print \$1}' snps.filtered.recode.imiss > remove_indv.txt
    # -------------------------------------------------------
    # 4) Remove individuals with high missingness from VCF
    # -------------------------------------------------------
    if [ -s remove_indv.txt ]; then
        vcftools --vcf snps.filtered.recode.vcf --remove remove_indv.txt --recode --out snps.filtered.rm_miss_indv
    else
        vcftools --vcf snps.filtered.recode.vcf --recode --out snps.filtered.rm_miss_indv
    fi
    # -------------------------------------------------------
    # 5) Remove sites with high missingness
    # -------------------------------------------------------
    vcftools --vcf snps.filtered.rm_miss_indv.recode.vcf --max-missing ${params.vcftools_site_max_miss} --recode --out snps.filtered.rm_miss_indv.site

    """
}
