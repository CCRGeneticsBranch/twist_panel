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

### 1. Point config to your data and output directories

The pipeline code, raw data, and output each live in separate directories.
Set absolute paths in `config/config.yaml`:

```yaml
raw_data_dir: "/path/to/raw/fastqs"   # your input fastq files
results_dir:  "/path/to/output"        # pipeline outputs (logs go in results_dir/logs/)
```

### 2. Prepare Input Data

Your fastq files stay in their existing location on the filesystem. Set
`raw_data_dir` to that directory in `config/config.yaml`.

### 3. Configure Pipeline

Edit `config/config.yaml`:

```yaml
samples:
  sample1:
    R1: "sample1_R1.fastq.gz"   # filename only - within raw_data_dir
    R2: "sample1_R2.fastq.gz"
  sample2:
    R1: "sample2_R1.fastq.gz"
    R2: "sample2_R2.fastq.gz"

raw_data_dir: "/path/to/raw/fastqs"
results_dir:  "/path/to/output"
reference_genome: "/fdb/bwa/indexes/hg38.fa"
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

### 4. Use HPC Reference Genome

The reference is pre-indexed on the HPC filesystem. Set `reference_genome` in `config/config.yaml` to that shared path (current default: `/fdb/bwa/indexes/hg38.fa`).

Optional check:

```bash
ls -lh /fdb/bwa/indexes/hg38.fa*
```

### 5. Run Pipeline

The pipeline is designed to run via the SLURM scheduler. Edit
`config/slurm_profile.yaml` to set your account and partition, then:

```bash
PIPELINE=/data/khanlab3/gb_omics/Erica/software/twist_panel

# Dry run — check what will execute without submitting any jobs
snakemake --snakefile $PIPELINE/Snakefile -n

# Submit all jobs to SLURM (standard usage)
snakemake --snakefile $PIPELINE/Snakefile \
    --profile $PIPELINE/config/slurm_profile.yaml

# Resume after failure (re-submits only incomplete jobs)
snakemake --snakefile $PIPELINE/Snakefile \
    --profile $PIPELINE/config/slurm_profile.yaml \
    --rerun-incomplete

# Run a single rule locally for testing (bypasses SLURM)
snakemake --snakefile $PIPELINE/Snakefile --use-envmodules fastqc_raw -j 1
```

### SLURM profile

Edit `config/slurm_profile.yaml` to match your cluster before the first run:

```yaml
jobs: 100                    # max simultaneous SLURM jobs
use-envmodules: true

default-resources:
  slurm_account: "your_account"   # ← change this
  slurm_partition: "norm"         # ← change this
  runtime: 240                    # minutes
  mem_mb: 8000
  cpus_per_task: 8
```

Per-rule overrides for resource-intensive steps are already set in
`config/slurm_profile.yaml` — adjust as needed.

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
- `reference_genome` points to the shared HPC FASTA path
- BWA index sidecar files exist at the same prefix
- `raw_data_dir` points to the directory containing your fastq files

```bash
# Check reference and index sidecars
ls -lh /fdb/bwa/indexes/hg38.fa*

# Check raw data
ls -lh /path/to/raw/fastqs/
```

### 3. Memory Issues

Increase `mem_mb` for the affected rule in `config/slurm_profile.yaml`:

```yaml
set-resources:
  merge_bam_alignment:
    mem_mb: 32000
```

Or raise the Java heap globally:

```yaml
memory_gb: 8   # in config/config.yaml
```

### 4. Module Loading Issues

Check available modules on your cluster:

```bash
module avail fastqc multiqc bwa samtools picard fgbio
```

Update the `modules` section in `config/config.yaml` with the exact names and versions shown.

## Tool Versions

The pipeline uses:
- **snakemake** ≥ 7.0
- **FastQC** 0.11.9
- **MultiQC** 1.12
- **Trimmomatic** (bundled)
- **BWA** 0.7.17
- **SAMtools** 1.15
- **umi_tools** 1.1.2

Update module versions in `config/config.yaml` under `modules` if needed.

## References

- Twist Bioscience: [DOC-001337 - Processing Sequencing Data Utilizing UMI](https://www.twistbioscience.com/sites/default/files/resources/2023-03/DOC-001337_TechNote-ProcessingSequencingDataUtilizingUMI-REV1-singles.pdf)
- umi_tools: [Documentation](https://umi-tools.readthedocs.io/)
- Snakemake: [Documentation](https://snakemake.readthedocs.io/)

## Support

For issues or questions, check:
1. `<results_dir>/logs/` directory for detailed error messages
2. The tool-specific documentation
3. The `summary_report.txt` for pipeline statistics

## License

This pipeline follows the same practices as the Twist Bioscience technical documentation.
