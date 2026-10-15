# Twist Bioscience Panel Preprocessing Pipeline

Snakemake pipeline for paired-end Illumina data from Twist Bioscience hybrid-capture panels with duplex UMIs. The default panel is a custom NCI_RhadoSeqFoc_TE panel. Pipeline calls duplex consensus reads by UMI with fgbio, aligns to the genome, and reports capture QC before and after consensus calling. Based on Twist TechNote DOC-001337 (*Processing Sequencing Data Utilizing UMI*, REV1).

## Pipeline steps

| Step | Rule | Tool | Description |
|---|---|---|---|
| 1 | `fastqc_raw` | FastQC | Raw read QC (R1 and R2) |
| 2 | `align_raw_fastq` | BWA-MEM, samtools | Align raw read pairs |
| 3 | `add_umi_tags_from_read_name` | Python ([tag_umi_from_read_name.py](workflow/scripts/tag_umi_from_read_name.py)), Picard | Parse UMIs from read names into `ZA`, `ZB` and `RX` tags |
| 4 | `collect_hs_metrics` | Picard `CollectHsMetrics` | Capture QC before consensus |
| 5 | `set_mate_information` | fgbio `SetMateInformation` | Fix mate information |
| 6 | `group_reads_by_umi` | fgbio `GroupReadsByUmi` | Group reads into families by mapped position and UMI |
| 7 | `call_duplex_consensus` | fgbio `CallDuplexConsensusReads` | Build duplex consensus reads |
| 8 | `realign_consensus` | Picard `SamToFastq`, BWA-MEM | Realign consensus reads |
| 9 | `merge_consensus_bam` | Picard `MergeBamAlignment` | Restore consensus tags on the aligned consensus reads |
| 10 | `collect_hs_metrics_post` | Picard `CollectHsMetrics` | Capture QC after consensus |
| 11 | `multiqc` | MultiQC | Combined QC report |

## Requirements

- Runs on NIH [Biowulf](https://hpc.nih.gov/) (SLURM and environment modules).
- Must have access to dependencies listed in [config/config.yaml](config/config.yaml).

## Installation

```bash
git clone https://github.com/CCRGeneticsBranch/twist_panel.git
```

## Quick start

The pipeline is designed to be run in a separate directory (/path/to/my_run).

```bash
# 1. Create a work folder with template config files
bash twist_panel/run_twist_panel --runmode=init --workdir=/path/to/my_run [--samplesheet=/path/to/samples.tsv]

# 2. Edit /path/to/my_run/config/config.yaml and config/samples.tsv (see Inputs)

# 3. Dry run (writes /path/to/my_run/dryrun.log)
bash twist_panel/run_twist_panel --runmode=dryrun --workdir=/path/to/my_run

# 4. Submit to SLURM
bash twist_panel/run_twist_panel --runmode=run --workdir=/path/to/my_run

# Unlock the work folder after an interrupted run
bash twist_panel/run_twist_panel --runmode=unlock --workdir=/path/to/my_run
```

`-h` shows all options; `-v` prints the pipeline and launcher versions. Per-rule SLURM resources are in `WORKDIR/config/cluster.yaml` (copied from [resources/cluster.yaml](resources/cluster.yaml) at `init`).

## Inputs

### Sample sheet (`config/samples.tsv`)

| Column | Description | Example |
|---|---|---|
| `Sample` | Sample name (must be unique) | `sample1` |
| `R1` | Read 1 FASTQ file name, relative to `raw_data_dir` | `sample1_R1.fastq.gz` |
| `R2` | Read 2 FASTQ file name, relative to `raw_data_dir` | `sample1_R2.fastq.gz` |

### Expected FASTQ header format

Standard Twist libraries follow the read structure "5M2S+T" (5 bp UMI + 2 bp skip + template). This pipeline instead starts from FASTQ files where UMIs have already been extracted into the read name:

```
@VH01088:45:AAJ2V5CM5:1:1101:60230:1000:GTCGT+CAGGA 1:N:0:AATGAGGTAA+GTATAACCTG
```

- `VH01088:45:AAJ2V5CM5:1:1101:60230:1000`: instrument and position
- `GTCGT+CAGGA`: UMI1 and UMI2
- `1:N:0:AATGAGGTAA+GTATAACCTG`: standard Illumina read flags and index sequences

The tagging step writes `ZA:Z:UMI1`, `ZB:Z:UMI2` and `RX:Z:UMI1-UMI2`, which are used by fgbio for grouping and duplex consensus calling.

### Configuration (`config/config.yaml`)

| Key | Default | Description (source) | Edit? |
|---|---|---|---|
| `raw_data_dir` | `/path/to/raw/fastqs` | Folder containing the FASTQ files | Yes |
| `reference_genome` | `/fdb/genomebrowser/fasta/hg38.fa` | Reference FASTA (UCSC hg38; requires `.fai` and `.dict` for Picard) | No |
| `reference_bwa` | `/fdb/bwa/indexes/hg38.fa` | BWA index prefix (hg38) | No |
| `bait_intervals`, `target_intervals` | `PIPELINE_HOME/resources/NCI_RhadoSeqFoc_TE.*.interval_list` | Bait/target interval lists for Picard `CollectHsMetrics`. The provided files were created from Twist probe and target BED files using Picard `BedToIntervalList`. Source BED provenance: **TODO**; redistribution permission: **TODO** | If panel differs |
| `umi_pattern` | `(?P<umi_1>[ACGTN]+)\\+(?P<umi_2>[ACGTN]+)` | Regex with named groups `umi_1` and `umi_2` matching UMIs in read names | If header format differs |
| `min_map_quality` | `10` | `GroupReadsByUmi --min-map-q` | As needed |
| `min_input_base_quality` | `30` | `CallDuplexConsensusReads --min-input-base-quality` | As needed |
| `min_reads` | `"2 1 1"` | `CallDuplexConsensusReads --min-reads` (duplex, strand 1, strand 2) | As needed |
| `samples`, `results_dir`, `tools` | `WORKDIR/...`, `PIPELINE_HOME/...` | Work folder and launcher/repository paths (filled by `init`) | No |

#### UMI pattern

By default, UMIs of any length separated by `+` are extracted. For strict fixed-length matching, for example 5 bp + 5 bp:

```yaml
umi_pattern: "(?P<umi_1>.{5})\\+(?P<umi_2>.{5})"
```

#### Consensus calling

`min_reads` is passed to `fgbio CallDuplexConsensusReads --min-reads`:

- `"2 1 1"` (default): duplex consensus required, at least one read from each strand
- `"1 0 0"`: more permissive; single-strand consensus reads are kept

## Outputs

All paths are relative to the work folder.

| Path | Description |
|---|---|
| `results/qc/fastqc/<sample>_R1_fastqc.html`, `_R2_fastqc.html` | FastQC reports |
| `results/bam/<sample>.merged.bam` / `.bai` | Aligned raw reads with UMI tags and read groups |
| `results/bam/<sample>.family_size_histogram.txt` | UMI family sizes |
| `results/bam/<sample>.consensus_final.bam` / `.bai` | Aligned duplex consensus reads with consensus tags |
| `results/qc/hsmetrics/<sample>.pre_consensus_hs_metrics.txt`, `.pre_consensus_per_target_coverage.txt` | Capture metrics before consensus |
| `results/qc/hsmetrics/<sample>.consensus_hs_metrics.txt`, `.consensus_per_target_coverage.txt` | Capture metrics after consensus |
| `results/qc/multiqc_report.html` | Combined QC report |
| `results/pipeline_version.txt` | Pipeline version used for the run |
| `runinfo.yaml`, `snakemake.log`, `logs/` | Run information, Snakemake log and per-rule SLURM logs |

## Tools and references

Modules are listed in [resources/tools.yaml](resources/tools.yaml). Twist's recommended versions are noted.

| Tool | Module | Used in | Reference |
|---|---|---|---|
| Snakemake | `snakemake/7.19.1` | workflow | Mölder et al. 2021, *F1000Research* [doi:10.12688/f1000research.29032.2](https://doi.org/10.12688/f1000research.29032.2) |
| Twist UMI processing | — | method | Twist Bioscience TechNote DOC-001337, *Processing Sequencing Data Utilizing UMI* |
| FastQC | `fastqc` | `fastqc_raw` | [Babraham Bioinformatics](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) |
| BWA-MEM | `bwa/0.7.17` | `align_raw_fastq`, `realign_consensus` | Li 2013, [arXiv:1303.3997](https://arxiv.org/abs/1303.3997) |
| SAMtools | `samtools` | BAM processing | Danecek et al. 2021, *GigaScience* [doi:10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008) |
| Picard | `picard` | read groups, HsMetrics, `SamToFastq`, `MergeBamAlignment` | [Broad Institute Picard](https://broadinstitute.github.io/picard/) |
| fgbio | `fgbio/2.0.2` | UMI grouping and duplex consensus | [fulcrumgenomics/fgbio](https://github.com/fulcrumgenomics/fgbio) |
| MultiQC | `multiqc` | `multiqc` | Ewels et al. 2016, *Bioinformatics* [doi:10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354) |

## License

Original pipeline code and documentation are licensed under MIT; see [LICENSE](LICENSE). The bundled panel interval lists retain any applicable third-party terms. The MIT license does not grant redistribution rights for external or bundled third-party resources.

## Versioning

The pipeline version is in [VERSION](VERSION) and changes are listed in [CHANGELOG.md](CHANGELOG.md). Each run records the pipeline version, git commit and launcher version in `runinfo.yaml` and `results/pipeline_version.txt`. The launch script is a copy of [biowulf-snakemake-launcher](https://github.com/CCRGeneticsBranch/biowulf-snakemake-launcher) (`bin/launcher.sh`).

## Contact

Erica Pehrsson, Genetics Branch, Center for Cancer Research, National Cancer Institute. Please report problems through GitHub issues.
