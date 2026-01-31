#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process MULTIQC {
  tag "multiqc"
  publishDir "${params.outdir}/multiqc", mode: 'copy'

  input:
  path qc_files

  output:
  path "multiqc_report.html"
  path "multiqc_data"

  script:
  """
  mkdir -p qc_inputs
  cp -f ${qc_files} qc_inputs/ || true
  multiqc qc_inputs -o .
  """
}
