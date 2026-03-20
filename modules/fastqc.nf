#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process FASTQC {
  tag "$sample"
  publishDir "${params.outdir}/fastqc", mode: 'copy'
  errorStrategy 'ignore' // Continue with other samples even if one fails

  input:
  tuple val(sample), path(r1), path(r2)

  output:
  tuple val(sample), path("*_fastqc.zip"), path("*_fastqc.html")

  script:
  """
  fastqc -t ${task.cpus} ${r1} ${r2}
  """
}
