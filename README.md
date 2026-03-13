# Twist Bioscience UMI Panel Preprocessing Pipeline

A comprehensive Snakemake pipeline for preprocessing paired-end sequencing data with UMIs embedded in fastq headers, following the Twist Bioscience technical note (DOC-001337).

## Overview

This pipeline processes raw sequencing reads through the following steps:

1. **UMI Extraction**: Extracts UMI sequences from fastq headers
2. **Quality Control**: FastQC and MultiQC for read quality assessment
3. **Adapter Trimming**: Removes Illumina adapters and low-quality bases using Trimmomatic
4. **Reference Alignment**: Maps reads to reference genome using BWA-MEM
5. **UMI Deduplication**: Removes PCR duplicates based on UMI information using umi_tools
6. **Summary Report**: Generates comprehensive pipeline statistics

## Quick Start

### 1. Installation

```bash
# Clone or navigate to the pipeline directory
cd /data/khanlab3/gb_omics/Erica/software/twist_panel

# (Optional) Create conda environment for snakemake
conda create -n snakemake snakemake -c conda-forge -c bioconda
conda activate snakemake
```

### 2. Prepare Input Data

Create the following directory structure:

```
twist_panel/
├── data/
│   └── raw/
│       ├── sample1_R1.fastq.gz
│       ├── sample1_R2.fastq.gz
│       ├── sample2_R1.fastq.gz
│       └── sample2_R2.fastq.gz
├── reference/
│   ├── hg38.fa          (FASTA reference genome)
│   └── hg38.fa.*        (BWA index files)
└── [rest of pipeline structure]
```

### 3. Configure Pipeline

Edit `config/config.yaml`:

```yaml
samples:
  sample1:
    R1: "sample1_R1.fastq.gz"
    R2: "sample1_R2.fastq.gz"
  sample2:
    R1: "sample2_R1.fastq.gz"
    R2: "sample2_R2.fastq.gz"

reference_genome: "reference/hg38.fa"
reference_index: "reference/hg38.fa.bwt"
```

Or use a CSV file:

```
# metadata/samples.csv
sample,R1,R2
sample1,sample1_R1.fastq.gz,sample1_R2.fastq.gz
sample2,sample2_R1.fastq.gz,sample2_R2.fastq.gz
```

Then in config.yaml:
```yaml
samples: "metadata/samples.csv"
```

### 4. Download Reference Genome

```bash
# For human genome (hg38)
cd reference
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Create BWA index
bwa index hg38.fa
```

### 5. Run Pipeline

```bash
# Dry run (see what will be executed)
snakemake -n

# Run with 8 cores
snakemake -j 8

# Run with all available cores
snakemake -j all

# Run specific rule (e.g., just QC)
snakemake fastqc_raw -j 8

# Run with conda environment creation
snakemake -j 8 --use-conda

# Run with conda and specific profile (e.g., slurm)
snakemake -j 100 --use-conda --profile slurm
```

### 6. View Results

```
results/
├── 01_qc/
│   ├── fastqc/           # FastQC reports
│   └── multiqc_report.html
├── 02_umi_extracted/
│   ├── sample1_R1.fastq.gz
│   ├── sample1_R2.fastq.gz
│   └── sample1_umi_stats.txt
├── 03_trimmed/
│   ├── sample1_R1.trimmed.fastq.gz
│   ├── sample1_R2.trimmed.fastq.gz
│   └── sample1.trimming_report.txt
├── 04_aligned/
│   ├── sample1.bam
│   └── sample1.bam.bai
├── 05_dedup/
│   ├── sample1.dedup.bam
│   ├── sample1.dedup.bam.bai
│   └── sample1_dedup_stats.txt
└── summary_report.txt
```

## Input Data Format

### Expected Fastq Header Format

The pipeline expects fastq files with UMIs embedded in the sequence header:

```
@VH01088:45:AAJ2V5CM5:1:1101:60230:1000:GTCGT+CAGGA 1:N:0:AATGAGGTAA+GTATAACCTG
```

Breaking this down:
- `VH01088:45:AAJ2V5CM5:1:1101:60230:1000` - Instrument and position info
- `GTCGT+CAGGA` - **UMI (this is extracted)**
- `1:N:0:AATGAGGTAA+GTATAACCTG` - Standard Illumina quality flags

### UMI Pattern

By default, the pipeline extracts 5bp+5bp UMIs:
```regex
(?P<umi_1>.{5})\+(?P<umi_2>.{5})
```

Modify in `config.yaml` if your UMI structure differs:
```yaml
umi_pattern: "(?P<umi_1>.{6})\\+(?P<umi_2>.{6})"  # For 6bp+6bp
```

## Pipeline Customization

### Modify Trimming Parameters

Edit `config.yaml`:

```yaml
trimming:
  steps:
    - "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True"
    - "LEADING:20"
    - "TRAILING:20"
    - "SLIDINGWINDOW:4:20"
    - "MINLEN:50"
```

### Change Computational Resources

```yaml
threads: 16          # Cores for parallel jobs
memory_gb: 32        # Memory allocation
```

### Enable Additional Analyses

```yaml
analyses:
  fastqc: true
  multiqc: true
  variant_calling: false    # Set to true for SNP calling
  coverage_analysis: false  # Set to true for coverage depth
```

### Adjust UMI Deduplication

```yaml
deduplication:
  method: "directional"  # directional, percentile, or cluster
  edit_distance: 1       # Max edit distance for UMI grouping
  per_base: false        # Per-base vs per-position UMI grouping
  subset: 0              # Use only first N reads (0 = all)
```

## Key Differences from Standard Twist Protocol

The standard Twist protocol assumes UMIs are in separate read files (e.g., R1 contains UMI, R2 contains sample sequence).

This pipeline is adapted for:
- **UMIs in headers**: Extracts from the fastq header instead of separate files
- **UMI prepending**: Attaches UMI to read sequence for downstream umi_tools compatibility
- **Integrated workflow**: All steps (extraction → QC → alignment → dedup) in one pipeline

## Output Description

### 02_umi_extracted/
- `*_R1.fastq.gz`, `*_R2.fastq.gz`: Fastqs with UMI prepended and header annotated
- `*_umi_stats.txt`: UMI extraction statistics

### 03_trimmed/
- `*.trimmed.fastq.gz`: Quality and adapter-trimmed reads
- `*.trimming_report.txt`: Trimmomatic statistics

### 04_aligned/
- `*.bam`: Aligned reads in coordinate-sorted BAM format
- `*.bam.bai`: BAM index files

### 05_dedup/
- `*.dedup.bam`: UMI-deduplicated reads
- `*_dedup_stats.txt`: UMI deduplication statistics

### summary_report.txt
- Comprehensive pipeline statistics and QC metrics

## Troubleshooting

### 1. UMI Extraction Failing

Check your UMI pattern in the fastq header. Verify with:
```bash
zcat data/raw/sample1_R1.fastq.gz | head -1
```

Update `umi_pattern` in `config.yaml` if needed.

### 2. Alignment Issues

Ensure:
- Reference genome is in FASTA format
- BWA index files exist (created by `bwa index`)
- File paths in config.yaml are correct

```bash
# Check reference
ls reference/
# Should show: hg38.fa, hg38.fa.amb, hg38.fa.ann, hg38.fa.bwt, etc.
```

### 3. Memory Issues

Increase memory allocation in `config.yaml`:
```yaml
memory_gb: 64  # Increase as needed
```

Or limit threads:
```bash
snakemake -j 4  # Use fewer parallel jobs
```

### 4. Conda Issues

Force environment update:
```bash
snakemake -j 8 --use-conda --conda-frontend mamba  # Faster
```

## Tool Versions

The pipeline uses:
- **snakemake** ≥ 7.0
- **FastQC** 0.11.9
- **MultiQC** 1.12
- **Trimmomatic** (bundled)
- **BWA** 0.7.17
- **SAMtools** 1.15
- **umi_tools** 1.1.2

Update in `envs/*.yaml` files if needed.

## References

- Twist Bioscience: [DOC-001337 - Processing Sequencing Data Utilizing UMI](https://www.twistbioscience.com/sites/default/files/resources/2023-03/DOC-001337_TechNote-ProcessingSequencingDataUtilizingUMI-REV1-singles.pdf)
- umi_tools: [Documentation](https://umi-tools.readthedocs.io/)
- Snakemake: [Documentation](https://snakemake.readthedocs.io/)

## Support

For issues or questions, check:
1. `logs/` directory for detailed error messages
2. The tool-specific documentation
3. The `summary_report.txt` for pipeline statistics

## License

This pipeline follows the same practices as the Twist Bioscience technical documentation.
