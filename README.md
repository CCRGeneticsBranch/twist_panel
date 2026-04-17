# Twist Bioscience UMI Panel Preprocessing Pipeline

A Snakemake pipeline for preprocessing paired-end sequencing data generated with the custom TWIST panel NCI_RhadoSeqFoc_TE 
Generally follows the Twist Bioscience technical note (DOC-001337 TechNote - Processing Sequencing Data Utilizing UMI), 
except UMIs are extracted into fastq headers during demultiplexing

## Overview

This pipeline processes raw sequencing reads through the following steps:

1. **FastQC**: Runs FastQC on raw R1 and R2 reads
2. **Primary Alignment**: Aligns raw fastqs with BWA-MEM
3. **UMI Tagging**: Copies UMI values from read names into ZA/ZB/RX BAM tags
4. **Pre-consensus QC**: Runs Picard CollectHsMetrics on merged BAMs
5. **Mate Fixing**: Runs fgbio SetMateInformation for paired-family grouping
6. **UMI Family Grouping**: Runs fgbio GroupReadsByUmi and writes family size histograms
7. **Consensus Calling**: Runs fgbio CallDuplexConsensusReads
8. **Consensus Realignment**: Converts to FASTQ and re-aligns consensus reads
9. **Consensus Merge**: Runs Picard MergeBamAlignment and restores consensus tags
10. **Post-consensus QC**: Runs Picard CollectHsMetrics on final consensus BAMs
11. **MultiQC**: Aggregates FastQC, both HsMetrics sets, and family size histograms

## Quick Start

Run the pipeline with [run_twist.sh](run_twist.sh) using `--runmode` and `--workdir`.

### 1. Initialize a run directory

This creates the work directory and copies the template config file and sample sheet.

```bash
bash run_twist.sh --runmode=init --workdir=/path/to/my_run
```

Then edit:
- `/path/to/my_run/config/config.yaml`
- `/path/to/my_run/config/samples.tsv`

### 2. Configure sample sheet and paths

Sample metadata must be TSV with columns `Sample`, `R1`, and `R2`.
`R1`/`R2` are filenames relative to `raw_data_dir`.

Example `samples.tsv`:

```tsv
Sample	R1	R2
sample1	sample1_R1.fastq.gz	sample1_R2.fastq.gz
sample2	sample2_R1.fastq.gz	sample2_R2.fastq.gz
```

Set raw data directory in `config.yaml`:

```yaml
raw_data_dir: "/path/to/raw/fastqs"
```

### 3. Dry run

Validate the rules without executing jobs.

```bash
bash run_twist.sh --runmode=dryrun --workdir=/path/to/my_run
```

### 4. Run on SLURM

Submit the pipeline to SLURM using the wrapper defaults and `config/biowulf.yaml`.

```bash
bash run_twist.sh --runmode=run --workdir=/path/to/my_run
```

### 5. Optional run modes

```bash
# Run interactively on an allocated node
bash run_twist.sh --runmode=runlocal --workdir=/path/to/my_run

# Unlock a stale Snakemake lock
bash run_twist.sh --runmode=unlock --workdir=/path/to/my_run

# Recreate workdir from scratch
bash run_twist.sh --runmode=reset --workdir=/path/to/my_run
```

### 6. Monitor run outputs

Key files written in the work directory:
- `snakemake.log`
- `snakemake.stats`
- `runinfo.yaml`
- `submit_script.sbatch`

### 7. View Results

```
results/
├── bam/
│   ├── sample1.merged.bam
│   ├── sample1.family_size_histogram.txt
│   └── sample1.consensus_final.bam
└── qc/
    ├── fastqc/
    ├── hsmetrics/
    │   ├── sample1.pre_consensus_hs_metrics.txt
    │   └── sample1.consensus_hs_metrics.txt
    ├── multiqc_report.html
    └── multiqc_data/
```

## Input Data Format

### Expected Fastq Header Format

Standard Twist libraries follow the format "5M2S+T" (5bp UMI + 2bp skip + template), 
and the protocol assumes UMI handling through fgbio read structures. 
But this pipeline starts from fastq files with UMIs embedded in the sequence header:

```
@VH01088:45:AAJ2V5CM5:1:1101:60230:1000:GTCGT+CAGGA 1:N:0:AATGAGGTAA+GTATAACCTG
```

Where:
- `VH01088:45:AAJ2V5CM5:1:1101:60230:1000` - Instrument and position info
- `GTCGT+CAGGA` - Extracted UMIs
- `1:N:0:AATGAGGTAA+GTATAACCTG` - Standard Illumina quality flags and barcodes

The custom tagging step parses UMI1 and UMI2 from read names and writes:
- ZA:Z:UMI1
- ZB:Z:UMI2
- RX:Z:UMI1-UMI2

These tags are then used by fgbio GroupReadsByUmi and downstream duplex consensus steps.

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

### Adjust UMI Grouping

Edit `config.yaml`:

```yaml
min_map_quality: 10
umi_edit_distance: 1
```

### Tune Consensus Calling

```yaml
error_rate_pre_umi: 45
min_reads: "2 1 1"
```

## Output Description

### bam/
- `*.merged.bam`, `*.merged.bai`: Primary aligned BAMs with UMI tags
- `*.family_size_histogram.txt`: Family size histograms from GroupReadsByUmi
- `*.consensus_final.bam`, `*.consensus_final.bai`: Final consensus alignments

### qc/fastqc/
- `*_fastqc.html`, `*_fastqc.zip`: FastQC reports for raw reads

### qc/hsmetrics/
- `*.pre_consensus_hs_metrics.txt`: Pre-consensus Picard HsMetrics
- `*.consensus_hs_metrics.txt`: Post-consensus Picard HsMetrics

### qc/
- `multiqc_report.html`: Aggregated report across FastQC, both HsMetrics sets, and family-size histograms
- `multiqc_data/`: Parsed data tables used to build the report

## Troubleshooting

### 1. UMI Tagging Failing

Check your UMI pattern in the fastq header. Verify with:
```bash
zcat data/raw/sample1_R1.fastq.gz | head -1
```

Update `umi_pattern` in `config.yaml` if needed.

### 2. Alignment Issues

Ensure:
- `reference_genome` and `reference_bwa` point to shared HPC FASTA paths
- BWA index sidecar files and dictionary exist at the same prefix
- `raw_data_dir` points to the directory containing your fastq files

```bash
# Check reference and index sidecars
ls -lh /fdb/bwa/indexes/hg38.fa*

# Check raw data
ls -lh /path/to/raw/fastqs/
```

### 3. Memory Issues

Increase `mem` for the affected rule in `config/biowulf.yaml`:

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
- **BWA** 0.7.17
- **SAMtools** 1.15
- **Picard** (module-provided)
- **fgbio** (module-provided)

Update module versions in `config/config.yaml` under `modules` if needed.

## References

- Twist Bioscience: [DOC-001337 - Processing Sequencing Data Utilizing UMI](https://www.twistbioscience.com/sites/default/files/resources/2023-03/DOC-001337_TechNote-ProcessingSequencingDataUtilizingUMI-REV1-singles.pdf)
- fgbio: [Documentation](https://fulcrumgenomics.github.io/fgbio/)

## Support

For issues or questions, check:
1. `<results_dir>/logs/` directory for detailed error messages
2. The tool-specific documentation
3. The MultiQC report and `multiqc_data/` tables for pipeline QC summaries

## License