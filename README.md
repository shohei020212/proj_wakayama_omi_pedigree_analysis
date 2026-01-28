# proj_wakayama_omi_pedigree_analysis

One-line description (e.g., A reproducible SNP calling and filtering pipeline using Nextflow and Docker).

## Overview
This repository provides a Nextflow (DSL2) pipeline to perform SNP analysis end-to-end (QC → mapping → variant calling → filtering → summary) in a reproducible manner. The pipeline is designed to run locally or on HPC/cloud by switching Nextflow configuration profiles. [web:9][web:26]

## Key features
- Workflow automation with Nextflow (restartable runs, scalable execution via profiles).
- Containerized execution with Docker to keep tool dependencies consistent across environments.
- Batch processing via a simple samplesheet (CSV).

---

## Requirements
- Nextflow
- Docker

> Note: Docker execution in Nextflow is controlled by settings such as `docker.enabled`.

---

## Quick start

### 1. Clone

```bash
$ git clone https://github.com/USER/REPO_NAME.git
$ cd REPO_NAME
```

###  2. Prepare inputs

Create samplesheet.csv:

```bash
sample,fastq_1,fastq_2
S1,data/S1_R1.fastq.gz,data/S1_R2.fastq.gz
S2,data/S2_R1.fastq.gz,data/S2_R2.fastq.gz
```


###  3. Run (Docker profile)

Create samplesheet.csv:

```bash
$ nextflow run main.nf -profile docker \
    --input samplesheet.csv \
    --fasta resources/reference.fa \
    --outdir results
```

## Pipeline steps

- Read QC (e.g., FastQC / fastp)
- Alignment to reference (e.g., BWA-MEM)
- Sorting / indexing (e.g., samtools)
- Variant calling (e.g., bcftools / GATK)
- Variant filtering (QUAL/DP/MAF, etc.)

## Summary statistics / reports

### Inputs
Minimum required:
```--input```: Path to a samplesheet CSV.
```--fasta```: Reference genome FASTA.
```--outdir```: Output directory.

Optional (examples):
```--gff```: Genome annotation (GFF/GTF).
```--bed```: Target regions.
```--known_sites```: Known variant sites for recalibration steps.

### Outputs
All outputs are written under --outdir (default: results/).
Example output layout:
results/qc/: QC reports
results/alignment/: BAM/CRAM (+ index)
results/variants/: Raw and filtered VCF/BCF
results/logs/: Process logs, pipeline logs
results/pipeline_info/: Execution metadata (versions/params, if implemented)

### Configuration profiles
Nextflow supports profiles (sets of configuration options) that can be selected at runtime using -profile.

Typical profiles (examples):
standard: Local execution (no containers)
docker: Docker-enabled execution
slurm / sge / pbs: HPC schedulers (if provided)