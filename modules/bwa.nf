#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process BWA {
  tag "$sample"
  publishDir "${params.outdir}/alignment", mode: 'symlink'

  cpus params.bwa_threads ? params.bwa_threads as int : 8
  memory params.bwa_mem ? params.bwa_mem : '16 GB'
  time params.bwa_time ? params.bwa_time : '8h'

  input:
  tuple val(sample), path(r1), path(r2)
  tuple path(ref_fasta), path(ref_index_files)

  output:
  tuple val(sample), path("${sample}.sam"), emit: sam

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
  """
}
