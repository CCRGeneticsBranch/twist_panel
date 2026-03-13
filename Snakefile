"""
Twist Bioscience UMI Panel Preprocessing Pipeline
Based on: DOC-001337 TechNote - Processing Sequencing Data Utilizing UMI (REV1)

Recommended tool versions (from the tech note):
  bwa        0.7.17
  samtools   1.15.1
  fgbio      1.5.1   (specifically required — newer versions change MergeBamAlignment)
  picard     2.27.1
  gatk       4.2.6.1

Standard Twist pipeline (from PDF)
-----------------------------------
1.  Picard   FastqToSam           raw fastqs → uBAM
2.  fgbio    ExtractUmisFromBam   inline UMIs (5M2S+T) → ZA / ZB / RX tags in uBAM
3.  Picard   SamToFastq           uBAM → interleaved fastq (UMI bases stripped)
4.  bwa mem -p                    interleaved fastq → coordinate-sorted aligned BAM
5.  Picard   MergeBamAlignment    merge aligned + UMI-tagged uBAM
6.  (opt)    CollectHsMetrics     pre-consensus coverage QC
7.  fgbio    GroupReadsByUmi      cluster read pairs by position + UMI (paired strategy)
8.  fgbio    CallDuplexConsensusReads  collapse families → duplex consensus uBAM
9.  Picard   SamToFastq           consensus uBAM → interleaved fastq
10. bwa mem -p                    realign consensus fastq
11. Picard   MergeBamAlignment    merge consensus aligned + consensus uBAM (retain tags)

Adaptation for data with UMIs already in the fastq header
-----------------------------------------------------------
Standard Twist library prep places UMI bases inline at the start of each read
(read structure 5M2S+T: 5 UMI bases, 2 skip bases, then template).
This data has UMIs already extracted into the read name header field:
  @VH01088:45:AAJ2V5CM5:1:1101:60230:1000:GTCGT+CAGGA 1:N:0:...
                                            ^^^^^  ^^^^^
                                            UMI1   UMI2

The custom pre-processing step (rule extract_umi_from_header) reads each
UMI pair from the header and prepends UMI1 to R1, UMI2 to R2, producing
fastqs with inline UMI prefixes compatible with the standard Twist pipeline.
Read structure is set to 5M+T (no skip bases, since they are not present in
the header-extracted UMIs). Update umi_read_structure in config.yaml if your
data does carry skip bases.
"""

import os
import pandas as pd
from pathlib import Path

configfile: "config/config.yaml"
# This workflow is configured for module-based software provisioning.
# Run Snakemake with: --use-envmodules

SAMPLES = config.get("samples", {})
if isinstance(SAMPLES, str):
    SAMPLES = pd.read_csv(SAMPLES, index_col=0).to_dict("index")

SAMPLE_NAMES = list(SAMPLES.keys()) if SAMPLES else []

RAW_DATA  = config.get("raw_data_dir",  "data/raw")
RESULTS   = config.get("results_dir",   "results")
LOGS      = config.get("logs_dir",      "logs")
REFERENCE = config.get("reference_genome", "")

THREADS = config.get("threads", 8)
MEMORY  = config.get("memory_gb", 4)   # per-tool heap (matches PDF: -Xmx4g)

# UMI parameters
UMI_PATTERN     = config.get("umi_pattern", r"(?P<umi_1>[ACGTN]+)\+(?P<umi_2>[ACGTN]+)")
# Read structure for fgbio ExtractUmisFromBam.
# Standard Twist inline: 5M2S+T  (5bp UMI + 2bp skip + template).
# For header-extracted UMIs prepended inline (no skip): 5M+T.
READ_STRUCTURE  = config.get("umi_read_structure", "5M+T")

# Pure module-based execution: set these in config.yaml to match your cluster.
MODULES = config.get("modules", {})
MOD_QC = MODULES.get("qc", ["fastqc", "multiqc"])
MOD_PREPROCESSING = MODULES.get("preprocessing", ["python"])
MOD_ALIGNMENT = MODULES.get("alignment", ["bwa", "samtools"])
MOD_PICARD = MODULES.get("picard", ["picard", "java"])
MOD_FGBIO = MODULES.get("fgbio", ["fgbio", "java"])
MOD_PICARD_ALIGN = MODULES.get(
    "picard_alignment",
    ["picard", "java", "bwa", "samtools"],
)


# ---------------------------------------------------------------------------
# Target rule
# ---------------------------------------------------------------------------

rule all:
    input:
        # FastQC / MultiQC on raw reads
        expand(f"{RESULTS}/01_qc/fastqc/{{sample}}_R1_fastqc.html", sample=SAMPLE_NAMES),
        expand(f"{RESULTS}/01_qc/fastqc/{{sample}}_R2_fastqc.html", sample=SAMPLE_NAMES),
        f"{RESULTS}/01_qc/multiqc_report.html",

        # Pre-processed fastqs (UMIs moved inline from header)
        expand(f"{RESULTS}/02_umi_extracted/{{sample}}_R1.fastq.gz", sample=SAMPLE_NAMES),
        expand(f"{RESULTS}/02_umi_extracted/{{sample}}_R2.fastq.gz", sample=SAMPLE_NAMES),

        # uBAM with UMI tags (ZA / ZB / RX)
        expand(f"{RESULTS}/04_ubam_umi_tagged/{{sample}}.umi_tagged.bam", sample=SAMPLE_NAMES),

        # Primary aligned BAM (after first MergeBamAlignment)
        expand(f"{RESULTS}/06_merged/{{sample}}.merged.bam", sample=SAMPLE_NAMES),
        expand(f"{RESULTS}/06_merged/{{sample}}.merged.bai", sample=SAMPLE_NAMES),

        # Grouped-by-UMI BAM
        expand(f"{RESULTS}/07_grouped/{{sample}}.grouped.bam", sample=SAMPLE_NAMES),

        # Final consensus BAM (aligned + UMI tags retained)
        expand(f"{RESULTS}/10_consensus_final/{{sample}}.consensus.final.bam",
               sample=SAMPLE_NAMES),
        expand(f"{RESULTS}/10_consensus_final/{{sample}}.consensus.final.bai",
               sample=SAMPLE_NAMES),

        f"{RESULTS}/summary_report.txt"


# ---------------------------------------------------------------------------
# Step 1 – FastQC / MultiQC on raw reads
# ---------------------------------------------------------------------------

rule fastqc_raw:
    """Run FastQC on raw fastq files."""
    input:
        lambda wc: f"{RAW_DATA}/{SAMPLES[wc.sample][wc.read]}"
    output:
        html = f"{RESULTS}/01_qc/fastqc/{{sample}}_{{read}}_fastqc.html",
        zip  = f"{RESULTS}/01_qc/fastqc/{{sample}}_{{read}}_fastqc.zip"
    params:
        outdir = f"{RESULTS}/01_qc/fastqc"
    log:
        f"{LOGS}/01_qc/fastqc_{{sample}}_{{read}}.log"
    envmodules:
        MOD_QC
    shell:
        "fastqc -o {params.outdir} -f fastq {input} 2>&1 | tee {log}"


rule multiqc:
    """Aggregate FastQC reports."""
    input:
        expand(f"{RESULTS}/01_qc/fastqc/{{sample}}_{{read}}_fastqc.zip",
               sample=SAMPLE_NAMES, read=["R1", "R2"])
    output:
        f"{RESULTS}/01_qc/multiqc_report.html"
    params:
        indir  = f"{RESULTS}/01_qc/fastqc",
        outdir = f"{RESULTS}/01_qc"
    log:
        f"{LOGS}/01_qc/multiqc.log"
    envmodules:
        MOD_QC
    shell:
        "multiqc {params.indir} -o {params.outdir} -f 2>&1 | tee {log}"


# ---------------------------------------------------------------------------
# Step 2 – Custom: move UMIs from header into reads inline
#
# Since UMIs are in the read name (not at the start of the read sequence)
# we prepend UMI1 → R1 and UMI2 → R2 so the standard fgbio ExtractUmisFromBam
# read-structure can strip them into ZA/ZB/RX tags downstream.
# ---------------------------------------------------------------------------

rule extract_umi_from_header:
    """
    Parse UMI1+UMI2 from the read-name and prepend each to its read:
      UMI1 → first 5 bases of R1  (read structure 5M+T)
      UMI2 → first 5 bases of R2  (read structure 5M+T)
    """
    input:
        r1 = lambda wc: f"{RAW_DATA}/{SAMPLES[wc.sample]['R1']}",
        r2 = lambda wc: f"{RAW_DATA}/{SAMPLES[wc.sample]['R2']}"
    output:
        r1    = f"{RESULTS}/02_umi_extracted/{{sample}}_R1.fastq.gz",
        r2    = f"{RESULTS}/02_umi_extracted/{{sample}}_R2.fastq.gz",
        stats = f"{RESULTS}/02_umi_extracted/{{sample}}_umi_stats.txt"
    params:
        script      = "scripts/extract_umi.py",
        umi_pattern = UMI_PATTERN
    log:
        f"{LOGS}/02_umi_extracted/{{sample}}.log"
    envmodules:
        MOD_PREPROCESSING
    shell:
        """
        python {params.script} \
            --r1 {input.r1} --r2 {input.r2} \
            --out_r1 {output.r1} --out_r2 {output.r2} \
            --stats {output.stats} \
            --umi_pattern '{params.umi_pattern}' \
            --split_umis \
            2>&1 | tee {log}
        """


# ---------------------------------------------------------------------------
# Step 3 – Picard FastqToSam  (PDF step 1)
# Convert paired fastqs to an unmapped BAM with sample/library metadata.
# ---------------------------------------------------------------------------

rule fastq_to_sam:
    """
    Picard FastqToSam: paired fastqs → unmapped BAM.

    PDF command:
      java -Xmx4g -jar picard.jar FastqToSam \\
        O={unaligned_bam} F1={r1} F2={r2} \\
        SM={sample} LB=Library1 PU=Unit1 PL=Illumina
    """
    input:
        r1 = f"{RESULTS}/02_umi_extracted/{{sample}}_R1.fastq.gz",
        r2 = f"{RESULTS}/02_umi_extracted/{{sample}}_R2.fastq.gz"
    output:
        ubam = temp(f"{RESULTS}/03_ubam/{{sample}}.unmapped.bam")
    params:
        mem     = MEMORY,
        sample  = lambda wc: wc.sample,
        library = lambda wc: config.get("library_name", wc.sample),
        platform = "Illumina"
    log:
        f"{LOGS}/03_ubam/{{sample}}.log"
    envmodules:
        MOD_PICARD
    shell:
        """
        picard --java-options "-Xmx{params.mem}g" FastqToSam \\
            O={output.ubam} \\
            F1={input.r1} \\
            F2={input.r2} \\
            SM={params.sample} \\
            LB={params.library} \\
            PU=Unit1 \\
            PL={params.platform} \\
            2>&1 | tee {log}
        """


# ---------------------------------------------------------------------------
# Step 4 – fgbio ExtractUmisFromBam  (PDF step 2)
# Strips the inline UMI prefix from each read and stores it in tags.
#   ZA = 5' UMI (from R1)
#   ZB = 3' UMI (from R2)
#   RX = combined "UMI1-UMI2"
# ---------------------------------------------------------------------------

rule extract_umis_from_bam:
    """
    fgbio ExtractUmisFromBam: inline UMI bases → ZA / ZB / RX BAM tags.

    PDF command:
      java -Xmx4g -jar fgbio.jar ExtractUmisFromBam \\
        --input={unaligned_bam} \\
        --output={unaligned_bam_umi_extracted} \\
        --read-structure=5M2S+T 5M2S+T \\
        --molecular-index-tags=ZA ZB \\
        --single-tag=RX

    Adaptation: read structure is 5M+T (no skip) because skip bases are not
    present in UMIs extracted from the read-name header.  Adjust
    umi_read_structure in config.yaml if your fastqs include skip bases.
    """
    input:
        ubam = f"{RESULTS}/03_ubam/{{sample}}.unmapped.bam"
    output:
        ubam_umi = f"{RESULTS}/04_ubam_umi_tagged/{{sample}}.umi_tagged.bam"
    params:
        rs    = READ_STRUCTURE
    log:
        f"{LOGS}/04_ubam_umi_tagged/{{sample}}.log"
    envmodules:
        MOD_FGBIO
    shell:
        """
        fgbio ExtractUmisFromBam \\
            --input={input.ubam} \\
            --output={output.ubam_umi} \\
            --read-structure={params.rs} {params.rs} \\
            --molecular-index-tags=ZA ZB \\
            --single-tag=RX \\
            2>&1 | tee {log}
        """


# ---------------------------------------------------------------------------
# Steps 5 + 6 – Picard SamToFastq  →  bwa mem -p  (PDF steps 3 + 4)
# Convert UMI-tagged uBAM → interleaved fastq, then align with bwa mem.
# ---------------------------------------------------------------------------

rule sam_to_fastq_and_align:
    """
    Picard SamToFastq (interleaved) piped directly into bwa mem -p.

    PDF commands:
      java -Xmx4g -jar picard.jar SamToFastq \\
        I={ubam_umi_extracted} F={fastq_umi_extracted} INTERLEAVE=true

      bwa mem -p -t 8 {ref} {fastq_umi_extracted} | \\
        samtools sort -@ 8 -o {aligned_bam}
    """
    input:
        ubam_umi = f"{RESULTS}/04_ubam_umi_tagged/{{sample}}.umi_tagged.bam",
        ref      = REFERENCE
    output:
        bam = temp(f"{RESULTS}/05_aligned/{{sample}}.aligned.bam")
    params:
        mem    = MEMORY
    threads: THREADS
    log:
        f"{LOGS}/05_aligned/{{sample}}.log"
    envmodules:
        MOD_PICARD_ALIGN
    shell:
        """
        picard --java-options "-Xmx{params.mem}g" SamToFastq \\
            I={input.ubam_umi} \\
            F=/dev/stdout \\
            INTERLEAVE=true \\
            QUIET=true \\
            2>>{log} | \\
        bwa mem -p -t {threads} {input.ref} /dev/stdin 2>>{log} | \\
        samtools sort -@ {threads} -o {output.bam} -
        """


# ---------------------------------------------------------------------------
# Step 7 – Picard MergeBamAlignment  (first merge, PDF step 5)
# Merges the aligned BAM with the UMI-tagged uBAM so that aligned records
# carry the ZA/ZB/RX tags and all Picard metadata.
# ---------------------------------------------------------------------------

rule merge_bam_alignment:
    """
    Picard MergeBamAlignment: restore UMI tags into the aligned BAM.

    PDF command:
      java -Xmx4g -jar picard.jar MergeBamAlignment \\
        UNMAPPED={unaligned_bam_umi_extracted} \\
        ALIGNED={aligned_bam_umi_extracted} \\
        O={aligned_tag_umi_bam} \\
        R={reference_fasta} \\
        CLIP_ADAPTERS=false \\
        VALIDATION_STRINGENCY=SILENT \\
        CREATE_INDEX=true \\
        EXPECTED_ORIENTATIONS=FR \\
        MAX_GAPS=-1 \\
        SO=coordinate \\
        ALIGNER_PROPER_PAIR_FLAGS=false
    """
    input:
        ubam    = f"{RESULTS}/04_ubam_umi_tagged/{{sample}}.umi_tagged.bam",
        aligned = f"{RESULTS}/05_aligned/{{sample}}.aligned.bam",
        ref     = REFERENCE
    output:
        bam = f"{RESULTS}/06_merged/{{sample}}.merged.bam",
        bai = f"{RESULTS}/06_merged/{{sample}}.merged.bai"
    params:
        mem    = MEMORY
    log:
        f"{LOGS}/06_merged/{{sample}}.log"
    envmodules:
        MOD_PICARD
    shell:
        """
        picard --java-options "-Xmx{params.mem}g" MergeBamAlignment \\
            UNMAPPED={input.ubam} \\
            ALIGNED={input.aligned} \\
            O={output.bam} \\
            R={input.ref} \\
            CLIP_ADAPTERS=false \\
            VALIDATION_STRINGENCY=SILENT \\
            CREATE_INDEX=true \\
            EXPECTED_ORIENTATIONS=FR \\
            MAX_GAPS=-1 \\
            SO=coordinate \\
            ALIGNER_PROPER_PAIR_FLAGS=false \\
            2>&1 | tee {log}
        """


# ---------------------------------------------------------------------------
# Step 8 (optional) – Picard CollectHsMetrics  (PDF step 6)
# Pre-consensus QC: fold-80, coverage, off-target rates.
# Requires a BED/interval file for your capture panel.
# ---------------------------------------------------------------------------

rule collect_hs_metrics:
    """
    Picard CollectHsMetrics: coverage QC on primary alignments.

    PDF command:
      java -Xmx4g -jar picard.jar CollectHsMetrics \\
        -I {aligned_tag_umi_bam} \\
        -O {metrics_file}	\\
        -R {reference_fasta} \\
        --BAIT_INTERVALS {probes_intervals} \\
        --TARGET_INTERVALS	target_file_UMI_demo_data_hg38.bed	\\		
        --PER_TARGET_COVERAGE	probe_file_UMI_demo_data_hg38.bed

    Only runs when a targets_intervals file is configured.
    """
    input:
        bam = f"{RESULTS}/06_merged/{{sample}}.merged.bam",
        ref = REFERENCE
    output:
        metrics = f"{RESULTS}/06_merged/{{sample}}.hs_metrics.txt"
    params:
        mem               = MEMORY,
        bait_intervals    = config.get("bait_intervals", ""),
        target_intervals  = config.get("target_intervals", ""),
        probe_file        = config.get("probe_file", "")
    log:
        f"{LOGS}/06_merged/{{sample}}.hs_metrics.log"
    envmodules:
        MOD_PICARD
    shell:
        """
        picard --java-options "-Xmx{params.mem}g" CollectHsMetrics \\
            -I {input.bam} \\
            -R {input.ref} \\
            -O {output.metrics} \\
            --BAIT_INTERVALS {params.bait_intervals} \\
            --TARGET_INTERVALS {params.target_intervals} \\
            --PER_TARGET_COVERAGE {params.probe_file} \\
            2>&1 | tee {log}
        """


# ---------------------------------------------------------------------------
# Step 9 – fgbio GroupReadsByUmi  (PDF step 7)
# Groups read pairs by mapped position + UMI.
# Must use --strategy=paired for CallDuplexConsensusReads.
# ---------------------------------------------------------------------------

rule group_reads_by_umi:
    """
    fgbio GroupReadsByUmi: cluster reads into UMI families.

    PDF command:
      java -Xmx4g -jar fgbio.jar GroupReadsByUmi \\
        --strategy=paired \\
        --input={aligned_tag_umi_bam} \\
        --output={grouped_by_umi_bam} \\
        --raw-tag=RX \\
        --min-map-q=10 \\
        --edits=1

    Note: --strategy=paired is required by CallDuplexConsensusReads.
    """
    input:
        bam = f"{RESULTS}/06_merged/{{sample}}.merged.bam",
        bai = f"{RESULTS}/06_merged/{{sample}}.merged.bai"
    output:
        bam  = f"{RESULTS}/07_grouped/{{sample}}.grouped.bam",
        hist = f"{RESULTS}/07_grouped/{{sample}}.family_size_histogram.txt"
    params:
        min_mapq  = config.get("min_map_quality", 10),
        edits     = config.get("umi_edit_distance", 1)
    log:
        f"{LOGS}/07_grouped/{{sample}}.log"
    envmodules:
        MOD_FGBIO
    shell:
        """
        fgbio GroupReadsByUmi \\
            --strategy=paired \\
            --input={input.bam} \\
            --output={output.bam} \\
            --family-size-histogram={output.hist} \\
            --raw-tag=RX \\
            --min-map-q={params.min_mapq} \\
            --edits={params.edits} \\
            2>&1 | tee {log}
        """


# ---------------------------------------------------------------------------
# Step 10 – fgbio CallDuplexConsensusReads  (PDF step 8)
# Assembles single-strand consensus reads, then combines
# the two strands into a duplex consensus read.
# ---------------------------------------------------------------------------

rule call_duplex_consensus:
    """
    fgbio CallDuplexConsensusReads: duplex error correction.

    PDF command:
      java -Xmx4g -jar fgbio.jar CallDuplexConsensusReads \\
        --input={grouped_by_umi_bam} \\
        --output={unaligned_consensus_bam} \\
        --error-rate-pre-umi=45 \\
        --error-rate-post-umi=30 \\
        --min-input-base-quality=30	\\
        --min-reads 2 1 1

    --min-reads 2 1 1:
      first integer  = duplex filter (≥2 total reads to form a consensus)
      second integer = min reads from the first strand
      third integer  = min reads from the second strand
    Twist recommend 2 1 1 to require evidence from both strands.
    """
    input:
        bam = f"{RESULTS}/07_grouped/{{sample}}.grouped.bam"
    output:
        ubam = f"{RESULTS}/08_consensus_ubam/{{sample}}.consensus.unmapped.bam"
    params:
        error_rate_pre    = config.get("error_rate_pre_umi", 45),
        min_reads         = config.get("min_reads", "2 1 1")
    log:
        f"{LOGS}/08_consensus_ubam/{{sample}}.log"
    envmodules:
        MOD_FGBIO
    shell:
        """
        fgbio CallDuplexConsensusReads \\
            --input={input.bam} \\
            --output={output.ubam} \\
            --error-rate-pre-umi={params.error_rate_pre} \\
            --min-reads {params.min_reads} \\
            2>&1 | tee {log}
        """


# ---------------------------------------------------------------------------
# Steps 11 + 12 – Picard SamToFastq  →  bwa mem -p  (PDF steps 9 re-alignment)
# Same pattern as steps 5+6, applied to consensus uBAM.
# ---------------------------------------------------------------------------

rule realign_consensus:
    """
    Picard SamToFastq + bwa mem -p: re-align consensus reads.

    PDF commands (same as steps 3+4):
      java -Xmx4g -jar picard.jar SamToFastq \\
        I={unaligned_consensus_bam} F={fastq_consensus} INTERLEAVE=true

      bwa mem -p -t 8 {ref} -o {aligned_bam_consensus} {fastq_consensus}
    """
    input:
        ubam = f"{RESULTS}/08_consensus_ubam/{{sample}}.consensus.unmapped.bam",
        ref  = REFERENCE
    output:
        bam  = temp(f"{RESULTS}/09_consensus_aligned/{{sample}}.consensus.aligned.bam")
    params:
        mem    = MEMORY
    threads: THREADS
    log:
        f"{LOGS}/09_consensus_aligned/{{sample}}.log"
    envmodules:
        MOD_PICARD_ALIGN
    shell:
        """
        picard --java-options "-Xmx{params.mem}g" SamToFastq \\
            I={input.ubam} \\
            F=/dev/stdout \\
            INTERLEAVE=true \\
            QUIET=true \\
            2>>{log} | \\
        bwa mem -p -t {threads} {input.ref} /dev/stdout 2>>{log} | \\
        samtools sort -@ {threads} -o {output.bam} -
        """


# ---------------------------------------------------------------------------
# Step 13 – Second Picard MergeBamAlignment  (PDF step 10)
# Merges the re-aligned consensus BAM with the consensus uBAM to restore
# duplex UMI tags (ZS, ZI, ZM, ZC, ZN, ...) needed for downstream tools.
# ---------------------------------------------------------------------------

rule merge_consensus_bam:
    """
    Picard MergeBamAlignment (second pass): retain UMI tag metadata in final BAM.

    PDF command:
      java -Xmx4g -jar picard.jar MergeBamAlignment \\
        VALIDATION_STRINGENCY=SILENT \\
        UNMAPPED={unaligned_consensus_bam} \\
        ALIGNED={aligned_bam_consensus} \\
        O={merged_consensus_bam} \\
        R={reference_fasta} \\
        CLIP_ADAPTERS=false \\
        CREATE_INDEX=true \\
        EXPECTED_ORIENTATIONS=FR \\
        MAX_GAPS=-1 \\
        SORT_ORDER=coordinate \\
        ALIGNER_PROPER_PAIR_FLAGS=false \\
        ATTRIBUTES_TO_RETAIN=X0 \\
        ATTRIBUTES_TO_RETAIN=ZS \\
        ATTRIBUTES_TO_RETAIN=ZI \\
        ATTRIBUTES_TO_RETAIN=ZM \\
        ATTRIBUTES_TO_RETAIN=ZC \\
        ATTRIBUTES_TO_RETAIN=ZN \\
        ATTRIBUTES_TO_RETAIN=ad \\  
        ATTRIBUTES_TO_RETAIN=bd \\ 
        ATTRIBUTES_TO_RETAIN=cd \\
        ATTRIBUTES_TO_RETAIN=ae \\
        ATTRIBUTES_TO_RETAIN=be \\
        ATTRIBUTES_TO_RETAIN=ce

        java -Xmx4g -jar /path/to/picard.jar AddOrReplaceReadGroups \\
            I={merged_no_read_group_bam} \\
            O={consensus_aligned_merged_bam} \\
            RGID={sample_id} \\
            RGLB={sample_lib} \\
            RGPL=Illumina \\
            RGSM={sample_name} \\
            RGPU=NA
    """
    input:
        ubam    = f"{RESULTS}/08_consensus_ubam/{{sample}}.consensus.unmapped.bam",
        aligned = f"{RESULTS}/09_consensus_aligned/{{sample}}.consensus.aligned.bam",
        ref     = REFERENCE
    output:
        bam = f"{RESULTS}/10_consensus_final/{{sample}}.consensus.final.bam",
        bai = f"{RESULTS}/10_consensus_final/{{sample}}.consensus.final.bai"
    params:
        mem    = MEMORY
    log:
        f"{LOGS}/10_consensus_final/{{sample}}.log"
    envmodules:
        MOD_PICARD
    shell:
        """
        picard --java-options "-Xmx{params.mem}g" MergeBamAlignment \\
            VALIDATION_STRINGENCY=SILENT \\
            UNMAPPED={input.ubam} \\
            ALIGNED={input.aligned} \\
            O={output.bam} \\
            R={input.ref} \\
            CLIP_ADAPTERS=false \\
            CREATE_INDEX=true \\
            EXPECTED_ORIENTATIONS=FR \\
            MAX_GAPS=-1 \\
            SORT_ORDER=coordinate \\
            ALIGNER_PROPER_PAIR_FLAGS=false \\
            ATTRIBUTES_TO_RETAIN=X0 \\
            ATTRIBUTES_TO_RETAIN=ZS \\
            ATTRIBUTES_TO_RETAIN=ZI \\
            ATTRIBUTES_TO_RETAIN=ZM \\
            ATTRIBUTES_TO_RETAIN=ZC \\
            ATTRIBUTES_TO_RETAIN=ZN \\
            2>&1 | tee {log}
        """

# ---------------------------------------------------------------------------
# Step 14 (optional) – Picard CollectHsMetrics  (PDF step 10a)
# Post-consensus QC: fold-80, coverage, off-target rates.
# Requires a BED/interval file for your capture panel.
# ---------------------------------------------------------------------------

rule collect_hs_metrics_post:
    """
    Picard CollectHsMetrics: coverage QC on primary alignments.

    PDF command:
      java -Xmx4g -jar picard.jar CollectHsMetrics \\
        -I {consensus_aligned_merged_bam} \\
        -O {metrics_file}	\\
        -R {reference_fasta} \\
        --BAIT_INTERVALS {probes_intervals} \\
        --TARGET_INTERVALS	target_file_UMI_demo_data_hg38.bed	\\		
        --PER_TARGET_COVERAGE	probe_file_UMI_demo_data_hg38.bed

    Only runs when a targets_intervals file is configured.
    """


# ---------------------------------------------------------------------------
# Summary report
# ---------------------------------------------------------------------------

rule generate_summary_report:
    """Collect per-sample statistics from all pipeline stages."""
    input:
        expand(f"{RESULTS}/02_umi_extracted/{{sample}}_umi_stats.txt",        sample=SAMPLE_NAMES),
        expand(f"{RESULTS}/07_grouped/{{sample}}.family_size_histogram.txt",   sample=SAMPLE_NAMES),
        expand(f"{RESULTS}/10_consensus_final/{{sample}}.consensus.final.bam", sample=SAMPLE_NAMES)
    output:
        report = f"{RESULTS}/summary_report.txt"
    log:
        f"{LOGS}/summary.log"
    envmodules:
        MOD_PREPROCESSING
    script:
        "scripts/generate_summary.py"
