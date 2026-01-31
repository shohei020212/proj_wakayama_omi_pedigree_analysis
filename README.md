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
- Java (version 17 or later)

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
sample,fastq1,fastq2
n1,data/n1.1.fq.gz,data/n1.2.fq.gz
n2,data/n2.1.fq.gz,data/n2.2.fq.gz
```


###  3. Run (Docker profile)

Create samplesheet.csv:

```bash
$ nextflow run main.nf \
    -profile docker \
    -params-file params.yaml \
    --input samplesheet.csv \
    --fasta resources/reference.fa \
    --outdir results \
    --trim_adapters resources/adapters.fa
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