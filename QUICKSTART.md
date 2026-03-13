# Quick Reference - Twist Panel Pipeline

## Directory Structure

```
twist_panel/
├── Snakefile                          # Main workflow definition
├── README.md                          # Full documentation
├── SETUP.md                           # Installation guide
├── ADVANCED.md                        # Advanced customization
│
├── config/
│   ├── config.yaml                   # Main configuration file
│   └── samples_template.csv           # Template for sample metadata
│
├── envs/
│   ├── preprocessing.yaml             # UMI extraction + utilities
│   ├── qc.yaml                        # FastQC, MultiQC
│   ├── alignment.yaml                 # BWA, samtools
│   └── dedup.yaml                     # umi_tools
│
├── scripts/
│   ├── extract_umi.py                 # UMI extraction script
│   └── generate_summary.py            # Summary report generator
│
├── data/
│   └── raw/                           # INPUT: Your fastq files here
│
├── reference/                         # INPUT: Reference genome files
│
└── results/                           # OUTPUT: Generated here
    ├── 01_qc/                         # FastQC/MultiQC reports
    ├── 02_umi_extracted/              # UMI-extracted fastqs
    ├── 03_trimmed/                    # Quality-trimmed fastqs
    ├── 04_aligned/                    # Aligned BAM files
    ├── 05_dedup/                      # Final deduplicated BAMs
    └── summary_report.txt             # Final statistics
```

## Most Common Commands

### Initial Setup

```bash
cd /data/khanlab3/gb_omics/Erica/software/twist_panel

# Install snakemake
conda create -n snakemake snakemake -c conda-forge

# Prepare directories
mkdir -p data/raw reference metadata

# Copy your fastqs
cp /path/to/*.fastq.gz data/raw/

# Activate environment
module load snakemake
```

### Download Reference (hg38)

```bash
cd reference
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
bwa index hg38.fa
cd ..
```

### Create Sample Metadata

```bash
# Create metadata/samples.csv with your sample names
cat > metadata/samples.csv << 'EOF'
sample,R1,R2
sample1,sample1_R1.fastq.gz,sample1_R2.fastq.gz
sample2,sample2_R1.fastq.gz,sample2_R2.fastq.gz
EOF
```

### Update Configuration

Edit `config/config.yaml`:
```yaml
samples: "metadata/samples.csv"
reference_genome: "reference/hg38.fa"
reference_index: "reference/hg38.fa.bwt"
threads: 8
```

### Run Full Pipeline

```bash
# Dry run (check for errors)
snakemake -n -j 1

# Run the pipeline
snakemake -j 8 --use-envmodules

# Run specific sample
snakemake extract_umi_from_header -j 8 --use-envmodules --config samples="metadata/sample1_only.csv"
```

## Output Files by Step

| Step | Output Files | Purpose |
|------|-------------|---------|
| 1 UMI Extract | `02_umi_extracted/*_R1.fastq.gz` | UMIs extracted from headers |
| 2 QC | `01_qc/multiqc_report.html` | FastQC quality report |
| 3 Trim | `03_trimmed/*_R1.trimmed.fastq.gz` | Quality-trimmed reads |
| 4 Align | `04_aligned/*.bam` | Aligned to reference genome |
| 5 Dedup | `05_dedup/*.dedup.bam` | PCR duplicates removed |
| 6 Report | `summary_report.txt` | Pipeline statistics |

## Key Files to Check

### Before Running
```bash
# Verify input fastqs
ls -lh data/raw/
zcat data/raw/sample1_R1.fastq.gz | head -2

# Verify reference
ls -lh reference/hg38.fa*

# Verify config
cat config/config.yaml
```

### After Running
```bash
# View summary
cat results/summary_report.txt

# Open QC report in browser
open results/01_qc/multiqc_report.html

# Check final BAMs
ls -lh results/05_dedup/*.dedup.bam

# Check alignment stats
samtools flagstat results/04_aligned/sample1.bam
```

## Troubleshooting Quick Fixes

### Problem: "Module conda not found"
**Fix**: `module load snakemake`

### Problem: "UMI extraction failed"
**Fix**: Check header format in first fastq line:
```bash
zcat data/raw/sample1_R1.fastq.gz | head -1
```
Update `umi_pattern` in config if format differs

### Problem: "Alignment failed - reference not found"
**Fix**:
```bash
ls -lh reference/
# Should show hg38.fa + multiple index files
```

### Problem: "Out of memory"
**Fix**: Reduce threads in config.yaml:
```yaml
threads: 4  # Was 8
```

## Configuration Template

**config/config.yaml minimum viable setup:**

```yaml
samples: "metadata/samples.csv"
raw_data_dir: "data/raw"
results_dir: "results"
logs_dir: "logs"
reference_genome: "reference/hg38.fa"
reference_index: "reference/hg38.fa.bwt"
threads: 8
memory_gb: 16
umi_pattern: "(?P<umi_1>.{5})\\+(?P<umi_2>.{5})"
```

## UMI Pattern Examples

Your fastq header format → regex pattern:

| Format | Pattern |
|--------|---------|
| `...GTCGT+CAGGA 1:N:0:...` (5+5) | `(?P<umi_1>.{5})\+(?P<umi_2>.{5})` |
| `...AGTCGT+CAGGA 1:N:0:...` (6+6) | `(?P<umi_1>.{6})\+(?P<umi_2>.{6})` |
| `...AGTCGTNN+CAGGA 1:N:0:...` (8+8) | `(?P<umi_1>.{8})\+(?P<umi_2>.{8})` |
| `...AGTCGT 1:N:0:...` (single) | `(?P<umi>[ACGT]{6})` |

## Performance Benchmarks

On 8-core system with 32GB RAM:

| Step | Time per Sample |
|------|-----------------|
| UMI extraction | 5-10 min |
| Trimming | 10-20 min |
| Alignment | 20-40 min |
| Deduplication | 5-15 min |
| **Total** | **~1-2 hours** |

## Important Notes

1. **UMI Quality**: Pipeline assumes UMIs are in header. Use `extract_umi.py` to validate extraction
2. **Reference**: Must be pre-indexed with `bwa index`
3. **Deduplication**: 50-80% duplicate removal is normal for PCR-amplified libraries
4. **Coverage**: Pipeline auto-calculates; no manual coverage computation needed
5. **BAM Sorting**: All output BAMs are coordinate-sorted and indexed

## Advanced Commands

```bash
# Visualize workflow DAG
snakemake --dag | dot -Tpng > dag.png

# Run with progress bar
snakemake -j 8 --use-envmodules --progress

# Rerun failed jobs only
snakemake -j 8 --use-envmodules --rerun-incomplete

# Use fewer cores for specific rule
snakemake -j 8 --use-envmodules -p --config threads=4

# Generate HTML report
snakemake -j 8 --use-envmodules --report results/workflow_report.html
```

## Getting Results Ready for Analysis

After pipeline completes:

```bash
# Copy final BAMs for downstream analysis
cp results/05_dedup/*.dedup.bam /path/to/analysis/
cp results/05_dedup/*.dedup.bam.bai /path/to/analysis/

# Generate per-sample statistics
for bam in results/05_dedup/*.dedup.bam; do
    echo "=== $(basename $bam) ==="
    samtools flagstat $bam
done

# Check mean coverage
for bam in results/05_dedup/*.dedup.bam; do
    cov=$(samtools depth $bam | awk '{sum+=$3; count++} END {print sum/count}')
    echo "$(basename $bam): ${cov}x"
done
```

## Support & Help

1. **Error Messages**: Check `logs/` directory for detailed error info
2. **Validation**: Run `snakemake -n` first (dry-run)
3. **Performance**: Adjust `threads` and `memory_gb` in config.yaml
4. **Custom Help**: See `ADVANCED.md` for advanced configurations

---

**For full documentation**: See README.md
**For setup steps**: See SETUP.md
**For customization**: See ADVANCED.md
