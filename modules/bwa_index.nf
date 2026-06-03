#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process BWA_INDEX {
  tag "bwa_index"
  publishDir "${params.outdir}/bwa_index", mode: 'copy'
  time = '24h'
  cpus   = 2

  input:
  path ref_fasta

  output:
  tuple path(ref_fasta), path("${ref_fasta}.*"), emit: ref_and_index

  script:
  """
  samtools faidx ${ref_fasta}
  bwa index ${ref_fasta}
  """
}
