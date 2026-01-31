#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process MULTIQC_FLAGSTAT {
  tag "multiqc_flagstat"
  publishDir "${params.outdir}/multiqc_flagstat", mode: 'copy'

  input:
  path flagstat_files
  
  output:
  path "multiqc_report.html"
  path "multiqc_data"

  script:
  """
  mkdir -p qc_inputs
  cp -f ${flagstat_files} qc_inputs/

  # QC report with MultiQC
  multiqc qc_inputs -o .
  """
}
