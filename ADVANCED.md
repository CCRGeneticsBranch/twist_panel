# Advanced Configuration & Tips

## Customizing UMI Extraction

### Different UMI Formats

Update `umi_pattern` in `config/config.yaml` based on your header format:

#### 5bp + 5bp (default)
```yaml
umi_pattern: "(?P<umi_1>.{5})\\+(?P<umi_2>.{5})"
```

#### 6bp + 6bp
```yaml
umi_pattern: "(?P<umi_1>.{6})\\+(?P<umi_2>.{6})"
```

#### 8bp + 8bp
```yaml
umi_pattern: "(?P<umi_1>.{8})\\+(?P<umi_2>.{8})"
```

#### Single barcode (no +)
```yaml
umi_pattern: "(?P<umi>.{8})"
```

#### UMI with specific sequence
```yaml
# If UMI starts with specific marker
umi_pattern: "BC:(?P<barcode>[ACGT]{8})"
```

### Test UMI Pattern

```bash
python3 << 'EOF'
import re
# Your header line
header = "@VH01088:45:AAJ2V5CM5:1:1101:60230:1000:GTCGT+CAGGA 1:N:0:AATGAGGTAA"

# Your pattern
pattern = r"(?P<umi_1>.{5})\+(?P<umi_2>.{5})"

match = re.search(pattern, header)
if match:
    print(f"✓ Pattern matches!")
    print(f"  UMI found: {match.group(0)}")
    print(f"  Groups: {match.groupdict()}")
else:
    print("✗ Pattern does not match")
    print("  Edit your regex pattern")
EOF
```

## Advanced Trimming Options

### Conservative Trimming (high quality)
```yaml
trimming:
  steps:
    - "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True"
    - "LEADING:25"
    - "TRAILING:25"
    - "SLIDINGWINDOW:4:25"
    - "MINLEN:75"
```

### Lenient Trimming (accept more reads)
```yaml
trimming:
  steps:
    - "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True"
    - "LEADING:15"
    - "TRAILING:15"
    - "SLIDINGWINDOW:4:15"
    - "MINLEN:35"
```

### No Quality Trimming (only adapters)
```yaml
trimming:
  steps:
    - "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True"
    - "MINLEN:35"
```

## Alternative Alignment Tools

### Using Bowtie2 Instead of BWA

Edit Snakefile to use bowtie2:

```python
rule align_to_reference:
    """Align with bowtie2."""
    input:
        r1 = f"{RESULTS}/03_trimmed/{{sample}}_R1.trimmed.fastq.gz",
        r2 = f"{RESULTS}/03_trimmed/{{sample}}_R2.trimmed.fastq.gz",
        ref = REFERENCE
    output:
        bam = f"{RESULTS}/04_aligned/{{sample}}.bam",
        bai = f"{RESULTS}/04_aligned/{{sample}}.bam.bai"
    threads: THREADS
    envmodules:
      "bowtie2/2.4.5"
      "samtools/1.15"
    shell:
        """
        bowtie2 -p {threads} -x {REFERENCE_INDEX} \
            -1 {input.r1} -2 {input.r2} | \
        samtools sort -@ {threads} -o {output.bam} - && \
        samtools index {output.bam}
        """
```

## Custom Deduplication Strategies

### Cluster-Based Deduplication (more conservative)
```yaml
deduplication:
  method: "cluster"
  edit_distance: 1
```

### Percentile-Based (for lower coverage)
```yaml
deduplication:
  method: "percentile"
  per_base: true
  subset: 100000  # Use first 100k reads for abundance estimation
```

### Custom UMI_Tools Parameters

Modify the umi_tools rule to include:

```python
rule umi_deduplication:
    # ... existing input/output ...
    shell:
        """
        umi_tools dedup -I {input.bam} \
            -S {output.dedup_bam} \
            --output-stats {params.temp_dedup} \
            --method cluster \
            --edit-distance-threshold 1 \
            --per-gene \
            --per-base \
            --log2stderr 2>{log} && \
        samtools index {output.dedup_bam}
        """
```

## Multi-Sample Patterns

### Process Samples in Batches

```yaml
# In config.yaml
batch_size: 5

# Then in Snakefile
BATCH_SAMPLES = [SAMPLE_NAMES[i:i+config['batch_size']] 
                 for i in range(0, len(SAMPLE_NAMES), config['batch_size'])]
```

### Process by Sample Group (e.g., different treatments)

```yaml
# config.yaml
sample_groups:
  treatment_a:
    - sample1
    - sample2
    - sample3
  treatment_b:
    - sample4
    - sample5
    - sample6
```

## Performance Optimization

### For Large Sample Numbers (>100 samples)

```bash
# Use job submission systems
snakemake -j 1000 --use-envmodules --cluster "qsub -l nodes=1:ppn=8 -l walltime=24:00:00"

# Or with SLURM
snakemake -j 1000 --use-envmodules --cluster "sbatch -c {threads} --mem=32GB"
```

### Parallel Rule Execution

```yaml
# config.yaml - Run multiple samples simultaneously
threads: 4
max_jobs: 8  # Run up to 8 samples in parallel
```

### Use Native Module Stacks

```bash
# Ask your HPC admins for a curated module set for this workflow
module avail fastqc multiqc bwa samtools picard fgbio
snakemake -j 8 --use-envmodules
```

## Quality Control Enhancements

### Generate Coverage Analysis

Add to Snakefile:

```python
rule coverage_analysis:
    """Calculate coverage depth statistics."""
    input:
        bam = f"{RESULTS}/05_dedup/{{sample}}.dedup.bam",
        bai = f"{RESULTS}/05_dedup/{{sample}}.dedup.bam.bai"
    output:
        coverage = f"{RESULTS}/06_qc/{{sample}}_coverage.txt"
    envmodules:
      "samtools/1.15"
    shell:
        """
        samtools depth {input.bam} > {output.coverage}
        """
```

### Add Per-Target Coverage

If you have a list of target regions:

```python
rule target_coverage:
    """Calculate coverage over specific targets."""
    input:
        bam = f"{RESULTS}/05_dedup/{{sample}}.dedup.bam",
        regions = "config/targets.bed"
    output:
        coverage = f"{RESULTS}/06_qc/{{sample}}_target_coverage.txt"
    envmodules:
      "bedtools/2.31"
    shell:
        """
        bedtools coverage -a {input.regions} -b {input.bam} > {output.coverage}
        """
```

## Variant Calling Addon

### Add Basic SNP Calling

```python
rule variant_calling:
    """Call variants using bcftools."""
    input:
        bam = f"{RESULTS}/05_dedup/{{sample}}.dedup.bam",
        ref = REFERENCE
    output:
        vcf = f"{RESULTS}/07_variants/{{sample}}.vcf.gz"
    envmodules:
      "bcftools/1.15"
      "samtools/1.15"
    shell:
        """
        samtools mpileup -u -f {input.ref} {input.bam} | \
        bcftools call -v -o {output.vcf} -
        """
```

## Checkpointing (Resume Failed Runs)

Snakemake automatically handles this, but you can control behavior:

```bash
# Resume interrupted run
snakemake -j 8 --use-envmodules --rerun-incomplete

# Force rerun of specific rule
snakemake -j 8 --use-envmodules --forcerun align_to_reference

# Continue from specific sample
snakemake -j 8 --use-envmodules --until extract_umi_from_header
```

## Debugging Failed Rules

### Increase Verbosity

```bash
# Show all commands
snakemake -j 1 --use-envmodules -rp

# Show output from failed commands
snakemake -j 1 --use-envmodules --keep-temp --keep-remote

# Print DAG visualizations
snakemake -j 1 --dag | dot -Tpng > dag.png
snakemake -j 1 --rulegraph | dot -Tpng > rulegraph.png
```

### Inspect Intermediate Files

```bash
# Snakemake keeps temp files with --keep-temp
snakemake -j 1 --keep-temp -rp

# Check what files were created
find results -type f | head -50
```

## Common Issues & Solutions

### Issue: Different UMI Counts in R1 vs R2

**Solution**: This is normal - each read pair is processed independently. The UMI is extracted from the same header for both reads.

### Issue: High Percentage of Duplicates Removed

**Solution**: This is expected with UMI-based approaches. 50-80% duplicate removal is typical for PCR-amplified libraries.

### Issue: Low Trimming Survival Rate

**Solution**: Your reads may have low quality. Try:
```yaml
trimming:
  steps:
    - "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True"
    - "LEADING:15"
    - "TRAILING:15"
    - "MINLEN:35"  # Reduce minimum length
```

### Issue: High Unmapped Read Rate

**Solution**: 
1. Verify reference genome is correct
2. Check read quality with FastQC report
3. Try less stringent trimming
4. Use `bwa mem -M` flag for split-read handling

## Useful One-Liners

```bash
# Check alignment statistics
samtools flagstat results/04_aligned/sample1.bam

# Check deduplication rate
samtools view -c results/04_aligned/sample1.bam
samtools view -c results/05_dedup/sample1.dedup.bam

# Calculate dedup percentage
echo "scale=2; $(samtools view -c results/05_dedup/sample1.dedup.bam) / $(samtools view -c results/04_aligned/sample1.bam) * 100" | bc

# Extract coverage statistics
samtools depth results/05_dedup/sample1.dedup.bam | awk '{sum+=$3; count++} END {print "Mean coverage: " sum/count}'

# Find unmapped reads
samtools view -f 4 results/04_aligned/sample1.bam | wc -l

# Check read length distribution
samtools view results/05_dedup/sample1.dedup.bam | cut -f10 | awk '{print length}' | sort | uniq -c | sort -rn | head
```
