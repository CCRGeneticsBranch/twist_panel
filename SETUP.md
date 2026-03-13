# Setup Instructions for Twist Panel UMI Pipeline

## Prerequisites

- Linux/macOS environment
- Conda (Miniconda or Anaconda)
- At least 50GB disk space (for reference genome and results)
- Internet connection for downloading reference genome

## Step-by-Step Setup

### 1. Clone or Download Pipeline

```bash
# Navigate to your software directory
cd /data/khanlab3/gb_omics/Erica/software/

# The pipeline directory is already set up at:
cd twist_panel
```

### 2. Install Snakemake

```bash
# Create a conda environment for snakemake
conda create -n snakemake snakemake=7.25 -c conda-forge -c bioconda

# Activate the environment
module load snakemake

# Verify installation
snakemake --version
```

### 3. Prepare Directory Structure

```bash
# Navigate to pipeline root
cd /data/khanlab3/gb_omics/Erica/software/twist_panel

# Create data and reference directories
mkdir -p data/raw
mkdir -p reference
mkdir -p metadata

# Directories for results and logs will be created automatically by Snakemake
```

### 4. Copy Your Input Data

```bash
# Copy your fastq files to data/raw/
# Example:
cp /path/to/your/fastqs/*.fastq.gz data/raw/

# List to verify
ls -lh data/raw/
```

### 5. Download and Prepare Reference Genome

#### Option A: Human Genome (hg38) - most common

```bash
cd reference

# Download reference genome
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# Extract
gunzip hg38.fa.gz

# Create BWA index (this takes ~2-3 hours for hg38)
module load snakemake
bwa index hg38.fa

# Verify index creation
ls -lh hg38.fa*
# Should show: hg38.fa.amb, hg38.fa.ann, hg38.fa.bwt, hg38.fa.pac, hg38.fa.sa

cd ..
```

#### Option B: Human Genome (GRCh38) - Ensembl version

```bash
cd reference

# Download from Ensembl
wget http://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Extract
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Create index
bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa

cd ..
```

#### Option C: Use Conda to get reference

```bash
# Use bioconda for pre-indexed genomes (faster)
conda install -c bioconda bwa-mem2

# Or download pre-indexed references (if available from your institution)
```

### 6. Create Sample Metadata File

Create `metadata/samples.csv`:

```csv
sample,R1,R2
sample1,sample1_R1.fastq.gz,sample1_R2.fastq.gz
sample2,sample2_R1.fastq.gz,sample2_R2.fastq.gz
sample3,sample3_R1.fastq.gz,sample3_R2.fastq.gz
```

Or use the template:

```bash
cp config/samples_template.csv metadata/samples.csv
# Edit metadata/samples.csv with your sample names
```

### 7. Update Configuration File

Edit `config/config.yaml` to point to your samples and reference:

```yaml
# Load samples from CSV
samples: "metadata/samples.csv"

# Or define inline:
# samples:
#   sample1:
#     R1: "sample1_R1.fastq.gz"
#     R2: "sample1_R2.fastq.gz"

# Update reference paths
reference_genome: "reference/hg38.fa"
reference_index: "reference/hg38.fa.bwt"  # Path to BWA index prefix

# Adjust compute resources as needed
threads: 8
memory_gb: 16
```

### 8. Verify Input Data Format

Check that your fastq headers contain UMIs in the expected format:

```bash
# Look at first header line
zcat data/raw/sample1_R1.fastq.gz | head -1

# Expected format:
# @VH01088:45:AAJ2V5CM5:1:1101:60230:1000:GTCGT+CAGGA 1:N:0:AATGAGGTAA+GTATAACCTG
#                                              ^^^^^^^^^
#                                              UMI portion
```

If your UMI format differs, update `umi_pattern` in `config/config.yaml`.

### 9. Test Configuration

```bash
# Activate snakemake environment
module load snakemake

# Run a dry-run to check for errors
snakemake -n -j 1

# Example output:
# Job counts:
#     count	jobs
#     1	all
#     1	extract_umi_from_header
#     ...
#     Total	XX
#
# This will execute total of XX jobs with cores: 1
```

### 10. Run Pipeline

```bash
# Ensure snakemake environment is active
module load snakemake

# Run with 8 cores (adjust based on your system)
snakemake -j 8 --use-envmodules

# Run with all available cores
snakemake -j all --use-envmodules

# Run specific sample only
snakemake -j 8 --use-envmodules --config samples=metadata/sample1_only.csv

# Run with verbose output
snakemake -j 8 --use-envmodules -rp
```

## Optional: Set Up Slurm/HPC Compatibility

For clusters with SLURM, create `config/slurm_profile.yaml`:

```yaml
jobs: 100
cores: 1024
account: "your_account"
partition: "normal"
default-resources:
  slurm_account: "your_account"
  slurm_partition: "normal"
  slurm_time: 1440
  slurm_mem_mb: 4000
```

Then run with:
```bash
snakemake -j 100 --use-envmodules --profile config/slurm_profile.yaml
```

## Troubleshooting

### Issue: "module: command not found"

**Solution**: Activate conda or ensure it's in your PATH
```bash
# If using Miniconda
source ~/miniconda3/bin/activate base
# or
export PATH="~/miniconda3/bin:$PATH"
```

### Issue: Reference genome download fails

**Solution**: Try alternative sources
```bash
# Use FTP instead of HTTP
ftp ftp.ensembl.org
# Navigate to: /pub/release-108/fasta/homo_sapiens/dna/

# Or use AWS S3
aws s3 cp s3://1000genomes/technical/reference/GRCh38/ reference/ --recursive
```

### Issue: "Could not install package"

**Solution**: Update conda
```bash
conda update conda
conda update mamba  # If using mamba backend
```

### Issue: Out of memory errors

**Solution**: Reduce threads in config.yaml
```yaml
threads: 4  # Instead of 8
memory_gb: 32
```

Or run with fewer parallel jobs:
```bash
snakemake -j 4 --use-envmodules
```

### Issue: UMI extraction fails

**Solution**: Check UMI pattern
```bash
# Inspect header format
zcat data/raw/sample1_R1.fastq.gz | sed -n '1p'

# Verify pattern in Python
python3 << 'EOF'
import re
header = "@VH01088:45:AAJ2V5CM5:1:1101:60230:1000:GTCGT+CAGGA 1:N:0:AATGAGGTAA+GTATAACCTG"
pattern = r"(?P<umi_1>.{5})\+(?P<umi_2>.{5})"
match = re.search(pattern, header)
print(f"Match: {match.group(0) if match else 'No match'}")
EOF
```

## Monitoring Pipeline Progress

```bash
# Check current status
snakemake -j 8 --use-envmodules --list-changes

# View DAG (dependency graph)
snakemake --dag | dot -Tpng > dag.png

# Real-time monitoring (in tmux/screen session)
watch -n 5 'ls -lh results/*/sample1* | tail -20'
```

## After Pipeline Completion

All results are in `results/` directory:

```bash
# View summary report
cat results/summary_report.txt

# Open MultiQC report in browser
open results/01_qc/multiqc_report.html

# List all final BAM files
ls results/05_dedup/*.dedup.bam
```

## Next Steps

Once deduplication is complete, you can:

1. **SNP/INDEL Calling** - Use GATK or bcftools:
   ```bash
   samtools mpileup -u -f reference/hg38.fa results/05_dedup/sample1.dedup.bam | bcftools call -v -o sample1.vcf.gz
   ```

2. **Coverage Analysis**:
   ```bash
   samtools depth results/05_dedup/sample1.dedup.bam | awk '{sum+=$3; count++} END {print sum/count}'
   ```

3. **IGV Viewer** - Visualize alignment:
   ```bash
   # BAM and BAI files are ready for IGV
   # They're in: results/05_dedup/*.dedup.bam*
   ```

## Performance Benchmarks

Typical run times on 8-core system with 32GB RAM:

- UMI extraction: 5-10 minutes per sample
- Adapter trimming: 10-20 minutes per sample
- Alignment: 20-40 minutes per sample
- UMI deduplication: 5-15 minutes per sample

**Total per sample**: ~1-2 hours

For multiple samples, batching is recommended.

## Getting Help

If you encounter issues:

1. Check the logs in `logs/` directory:
   ```bash
   tail -100 logs/04_aligned/sample1.log
   ```

2. Run with verbose output:
   ```bash
   snakemake -j 1 --use-envmodules -rp 2>&1 | head -100
   ```

3. Verify tool installations:
   ```bash
   module load snakemake
   fastqc --version
   bwa
   samtools --version
   umi_tools --version
   ```
