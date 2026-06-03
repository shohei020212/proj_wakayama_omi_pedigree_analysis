#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process REPAIR_FASTQ {
    tag "$sample"
    time   = { 1.m * task.attempt }
    cpus   = 2
    errorStrategy 'ignore'

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    tuple val(sample), path("${sample}_fixed.1.fastq.gz"), path("${sample}_fixed.2.fastq.gz"), emit: repaired_reads

    script:
    """
    # Download BBMap repair script
    wget https://github.com/BioInfoTools/BBMap/blob/master/sh/repair.sh
    chmod +x repair.sh
    
    # Fix files of paired reads that became disordered using BBMap
    repair.sh \
        in1=${r1} \
        in2=${r2} \
        out1=${sample}_fixed.1.fastq.gz \
        out2=${sample}_fixed.2.fastq.gz \
        outs=${sample}_singletons.fastq.gz \
        repair

    rm repair.sh
    """
}