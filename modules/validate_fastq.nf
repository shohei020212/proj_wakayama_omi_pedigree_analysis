#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process VALIDATE_FASTQ {
    tag "$sample"
    errorStrategy 'ignore' // Continue with other samples even if one fails validation

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    tuple val(sample), path(r1), path(r2), emit: validated_reads

    script:
    """
    echo "Validating $sample..."

    # 1. Get line counts for R1 and R2
    r1_lines=\$(zcat ${r1} | wc -l)
    r2_lines=\$(zcat ${r2} | wc -l)

    # 2. If R1 is empty, report an error and skip this sample
    if [ "\$r1_lines" -eq 0 ]; then
        echo "Error: $sample R1 is empty." >&2
        exit 1
    fi

    # 3. R1/R2 line count consistency check (should be the same for paired-end data)
    if [ "\$r1_lines" -ne "\$r2_lines" ]; then
        echo "Error: $sample line count mismatch (R1: \$r1_lines, R2: \$r2_lines)." >&2
        exit 1
    fi
    
    echo "$sample passed validation."
    """
}