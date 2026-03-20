#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process TRIMMOMATIC {
  tag "$sample"
  publishDir "${params.outdir}/trim", mode: 'copy'
  errorStrategy 'ignore' // Continue with other samples even if one fails

  input:
  tuple val(sample), path(r1), path(r2)
  path adapters

  output:

  tuple val(sample),
        path("${sample}.R1.trimmed.fq.gz"),
        path("${sample}.R2.trimmed.fq.gz"),
        emit: reads

  script:
  """
  trimmomatic PE -threads ${task.cpus} \
    ${r1} ${r2} \
    ${sample}.R1.trimmed.fq.gz \
    ${sample}.R1.unpaired.fq.gz \
    ${sample}.R2.trimmed.fq.gz \
    ${sample}.R2.unpaired.fq.gz \
    ILLUMINACLIP:${adapters}:2:30:10 \
    LEADING:${params.trim_leading} \
    TRAILING:${params.trim_trailing} \
    SLIDINGWINDOW:${params.trim_slidingwindow} \
    MINLEN:${params.trim_minlen} \
    ${params.trim_extra}
  """
}
