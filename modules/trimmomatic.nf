#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process TRIMMOMATIC {
  tag "$sample"
  publishDir "${params.outdir}/trim", mode: 'copy'

  cpus params.trim_threads as int

  input:
  tuple val(sample), path(r1), path(r2)
  path adapters

  output:

  tuple val(sample),
        path("${sample}.R1.trimmed.fq.gz"),
        path("${sample}.R2.trimmed.fq.gz"),
        emit: reads

  script:
  // Use only when using ILLUMINACLIP
  def adapterPart = adapters ? "ILLUMINACLIP:${adapters}:2:30:10" : ""
  """
  trimmomatic PE -threads ${task.cpus} \
    ${r1} ${r2} \
    ${sample}.R1.trimmed.fq.gz \
    ${sample}.R1.unpaired.fq.gz \
    ${sample}.R2.trimmed.fq.gz \
    ${sample}.R2.unpaired.fq.gz \
    ${adapterPart} \
    LEADING:20 \
    TRAILING:20 \
    SLIDINGWINDOW:4:20 \
    MINLEN:100 \
    ${params.trim_extra}
  """
}
