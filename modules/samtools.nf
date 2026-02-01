#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process SAMTOOLS_VIEW_SORT {
  tag "$sample"
  publishDir "${params.outdir}/alignment", mode: 'symlink'

  input:
  tuple val(sample), path(sam)
  
  output:
  tuple val(sample), path("${sample}_sorted.bam"), emit: bam
  tuple val(sample), path("${sample}_sorted.bam.bai"), emit: bai
  tuple val(sample), path("${sample}_filtered_flagstats.txt"), emit: flagstats_filtered
  tuple val(sample), path("${sample}_unfiltered_flagstats.txt"), emit: flagstats_unfiltered

  script:
  """
  # run samtools to convert SAM to filtered BAM + sort
  samtools view \
    -@ ${task.cpus} \
    -q ${params.samtools_Q ?: 0} \
    -F ${params.samtools_F ?: 0} \
    -bS \
    ${sam} \
    | samtools sort \
    -@ ${task.cpus} \
    -o ${sample}_sorted.bam

  # index the sorted BAM
  samtools index -@ ${task.cpus} ${sample}_sorted.bam

  # generate flagstats
  samtools flagstat -@ ${task.cpus} ${sample}_sorted.bam > ${sample}_filtered_flagstats.txt
  samtools flagstat -@ ${task.cpus} ${sample}.sam > ${sample}_unfiltered_flagstats.txt
  """
}
