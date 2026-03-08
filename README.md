# proj_wakayama_omi_pedigree_analysis

A reproducible SNP calling, filtering and pedigree recontruction pipeline using Nextflow and Docker.

##　🌏 Overview
This repository provides a Nextflow (DSL2) pipeline to perform SNP analysis end-to-end (QC → mapping → variant calling → filtering → QC → pedigree reconstruction → summary) in a reproducible manner. The pipeline is designed to run locally by switching Nextflow configuration profiles.

##　🔑 Key features
- Workflow automation with Nextflow (restartable runs, scalable execution via profiles).
- Containerized execution with Docker to keep tool dependencies consistent across environments.
- Batch processing via a simple samplesheet (CSV).

---

## 🛠 Prerequisites
- Nextflow
- Docker
- Java (version 17 or later)

> Note: Docker execution in Nextflow is controlled by settings such as `docker.enabled`.

---

## 🚀 Quick start

### 1. Clone

```bash
git clone https://github.com:shohei020212/proj_wakayama_omi_pedigree_analysis.git
cd proj_wakayama_omi_pedigree_analysis
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
nextflow run main.nf \
    -profile docker \
    -params-file params.yaml \
    --input samplesheet.csv \
    --fasta resources/reference.fa \
    --outdir results \
    --trim_adapters resources/adapters.fa \
    --regions resources/target_regions.bed \
    --known_sites resources/target_snplist.tsv \
    --lh_file resources/LH_Data.txt
```

## ⚙️ Pipeline steps

- FASTQ validation (Ubuntu)
- Read QC (FastQC / MultiQC)
- Adapter trimming (Trimmomatic)
- Alignment to reference (BWA-MEM)
- Sorting / indexing (SAMtools)
- Variant calling (bcftools or GATK)
- Variant filtering (VCFtools)

## 📊 Summary statistics / reports

### Inputs
Minimum required:  
```--input```: Path to a samplesheet CSV.  
```--fasta```: Reference genome FASTA.  
```--outdir```: Output directory.  
```--trim_adapters```: Adapter fasta file.   

Optional:  
```--caller```: Variant caller (defauls: bcftools, option: gatk and freebayes).  
```--regions```: Target region list (BED format: chr/start/end).   
```--known_sites```: Known variant sites (TSV; chr/pos).   
```--lh_file```: Life-history data (ID, Sex, BirthYear, BY.min, BY.max, Year.last).  

### Outputs
All outputs are written under `--outdir` (default: `results/`).  
Example output layout:  
`results/fastqc/`: QC reports  
`results/trim/`: trimmed FASTQ  
`results/alignment/`: SAM/BAM (+ index)  
`results/multiqc/`: FASTQC summary  
`results/multiqc_flagstat`: Alignment QC summary  
`results/bcftools(gatk or freebayes)`: Raw variants (VCF)  
`results/vcftools`: Filtered variants (VCF)  
`results/genotype_qc`: Genotype quality check based on allele counts (PNG)  
`results/sequoia`: Results of pedigree reconstruction  
`results/pipeline_info/`: resource-usage reports  

### Configuration profiles
Nextflow supports profiles (sets of configuration options) that can be selected at runtime using -profile.

Typical profiles (examples):
standard: Local execution (no containers)
docker: Docker-enabled execution

### Resource-usage reports
The pipeline automatically generates resource-usage reports after each run. These files are written to `results/pipeline_info/` and include a per-task execution trace (`trace.tsv`) as well as HTML summaries (`report.html` and `timeline.html`) that help you review CPU, memory, and runtime usage for each process.


