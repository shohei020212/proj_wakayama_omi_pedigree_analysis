#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process BWA_SAMTOOLS {
  tag "$sample"
  publishDir "${params.outdir}/alignment", mode: 'copy'

  input:
  tuple val(sample), path(r1), path(r2)
  tuple path(ref_fasta), path(ref_index_files)

  output:
  tuple val(sample), path("${sample}_sorted.bam"), emit: bam
  tuple val(sample), path("${sample}_sorted.bam.bai"), emit: bai
  tuple val(sample), path("${sample}_filtered_flagstats.txt"), emit: flagstats_filtered
  tuple val(sample), path("${sample}_unfiltered_flagstats.txt"), emit: flagstats_unfiltered


  script:
  """
  # Safety settings to make it "fail fast"
  set -euo pipefail

  # run bwa mem
  bwa mem \
    -t ${task.cpus} \
    -M \
    ${ref_fasta} \
    ${r1} \
    ${r2} \
    > ${sample}.sam

  # add read group information to BAM
  samtools addreplacerg \
    -r '@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA\\tLB:lib1' \
    -o ${sample}.rg.bam \
    ${sample}.sam
  
  # run samtools to convert SAM to filtered BAM + sort
  samtools view \
    -@ ${task.cpus} \
    -q ${params.samtools_Q ?: 0} \
    -F ${params.samtools_F ?: 0} \
    -bS \
    ${sample}.rg.bam \
    | samtools sort \
    -@ ${task.cpus} \
    -o ${sample}_sorted.bam

  # index the sorted BAM
  samtools index -@ ${task.cpus} ${sample}_sorted.bam

  # generate flagstats
  samtools flagstat -@ ${task.cpus} ${sample}_sorted.bam > ${sample}_filtered_flagstats.txt
  samtools flagstat -@ ${task.cpus} ${sample}.sam > ${sample}_unfiltered_flagstats.txt

  # cleanup intermediate files
  rm ${sample}.sam
  rm ${sample}.rg.bam
  """
}
