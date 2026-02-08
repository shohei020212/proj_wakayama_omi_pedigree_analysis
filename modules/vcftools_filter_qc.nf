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
    path "snps_filtered_${params.caller}.recode.vcf"
    path "snps_filtered_${params.caller}.rm_miss_indv.recode.vcf"
    path "snps_filtered_${params.caller}.rm_miss_indv.site.recode.vcf", emit: filtered_vcf

    // QC outputs
    path "snps_filtered_${params.caller}.recode.frq"
    path "snps_filtered_${params.caller}.recode.idepth"
    path "snps_filtered_${params.caller}.recode.ldepth.mean"
    path "snps_filtered_${params.caller}.recode.lqual"
    path "snps_filtered_${params.caller}.recode.imiss"
    path "snps_filtered_${params.caller}.recode.lmiss"
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
    --max-meanDP  ${params.vcftools_max_meanDP} \
    --recode \
    --out snps_filtered_${params.caller}
    # -------------------------------------------------------
    # 2) QC metrics
    # -------------------------------------------------------
    vcftools \
    --vcf snps_filtered_${params.caller}.recode.vcf \
    --freq2 \
    --out snps_filtered_${params.caller}.recode

    vcftools \
    --vcf snps_filtered_${params.caller}.recode.vcf \
    --depth \
    --out snps_filtered_${params.caller}.recode

    vcftools \
    --vcf snps_filtered_${params.caller}.recode.vcf \
    --site-mean-depth \
    --out snps_filtered_${params.caller}.recode

    vcftools \
    --vcf snps_filtered_${params.caller}.recode.vcf \
    --site-quality \
    --out snps_filtered_${params.caller}.recode

    vcftools \
    --vcf snps_filtered_${params.caller}.recode.vcf \
    --missing-indv \
    --out snps_filtered_${params.caller}.recode

    vcftools \
    --vcf snps_filtered_${params.caller}.recode.vcf \
    --missing-site \
    --out snps_filtered_${params.caller}.recode
    # -------------------------------------------------------
    # 3) Create a list of removed individuals due to missingness
    # -------------------------------------------------------
    awk '\$5 > ${params.vcftools_ind_max_miss} {print \$1}' snps_filtered_${params.caller}.recode.imiss > remove_indv.txt
    # -------------------------------------------------------
    # 4) Remove individuals with high missingness from VCF
    # -------------------------------------------------------
    if [ -s remove_indv.txt ]; then
        vcftools --vcf snps_filtered_${params.caller}.recode.vcf --remove remove_indv.txt --recode --out snps_filtered_${params.caller}.rm_miss_indv
    else
        vcftools --vcf snps_filtered_${params.caller}.recode.vcf --recode --out snps_filtered_${params.caller}.rm_miss_indv
    fi
    # -------------------------------------------------------
    # 5) Remove sites with high missingness
    # -------------------------------------------------------
    vcftools --vcf snps_filtered_${params.caller}.rm_miss_indv.recode.vcf --max-missing ${params.vcftools_site_max_miss} --recode --out snps_filtered_${params.caller}.rm_miss_indv.site

    """
}
