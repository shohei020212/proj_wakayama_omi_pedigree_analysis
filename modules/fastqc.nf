#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process FASTQC {
  tag "$sample"
  publishDir "${params.outdir}/fastqc", mode: 'copy'
  time   = { 1.m * task.attempt }
  cpus   = 4

  input:
  tuple val(sample), path(r1), path(r2)

  output:
  tuple val(sample), path("*_fastqc.zip"), path("*_fastqc.html")

  script:
  """
  fastqc -t ${task.cpus} ${r1} ${r2}
  """
}
