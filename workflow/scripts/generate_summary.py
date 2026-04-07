#!/usr/bin/env python3
"""
Generate a summary report for the Twist UMI panel pipeline (DOC-001337).
Collects statistics from: UMI extraction, UMI family size histograms,
and final consensus BAMs.
"""

import os
import glob
from pathlib import Path


def read_umi_stats(stats_file):
    stats = {"total_reads": 0, "umis_extracted": 0, "success_rate": 0.0}
    try:
        with open(stats_file) as f:
            for line in f:
                if "Total reads processed:" in line:
                    stats["total_reads"] = int(line.split(":")[1].strip())
                elif "UMIs successfully extracted:" in line:
                    stats["umis_extracted"] = int(line.split(":")[1].strip())
                elif "Extraction success rate:" in line:
                    stats["success_rate"] = float(
                        line.split(":")[1].strip().rstrip("%")
                    )
    except Exception as e:
        print(f"Warning: could not parse {stats_file}: {e}")
    return stats


def read_family_histogram(hist_file):
    """
    Parse fgbio GroupReadsByUmi family size histogram.
    Returns dict with median family size and counts.
    """
    stats = {"total_families": 0, "median_family_size": "N/A"}
    try:
        sizes = []
        with open(hist_file) as f:
            next(f)  # skip header
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    size, count = int(parts[0]), int(parts[1])
                    sizes.extend([size] * count)
        if sizes:
            stats["total_families"] = len(sizes)
            sizes.sort()
            mid = len(sizes) // 2
            stats["median_family_size"] = (
                sizes[mid] if len(sizes) % 2 != 0
                else (sizes[mid - 1] + sizes[mid]) / 2
            )
    except Exception as e:
        print(f"Warning: could not parse {hist_file}: {e}")
    return stats


def generate_report(output_file):
    results_dir = "results"

    umi_stats_files = glob.glob(f"{results_dir}/02_umi_extracted/*_umi_stats.txt")
    samples = sorted(
        Path(f).stem.replace("_umi_stats", "") for f in umi_stats_files
    )

    report = []
    report.append("=" * 80)
    report.append("TWIST PANEL UMI PREPROCESSING PIPELINE — SUMMARY REPORT")
    report.append("Based on: DOC-001337 TechNote (Picard + fgbio 1.5.1 workflow)")
    report.append("=" * 80)
    report.append("")

    if not samples:
        report.append("WARNING: No sample stats found in results directory.")
    else:
        hdr = f"{'Sample':<30} {'Total Reads':<18} {'UMI Extr%':<12} {'Families':<14} {'Median Family'}"
        report.append(hdr)
        report.append("-" * 80)

        for sample in samples:
            umi  = read_umi_stats(f"{results_dir}/02_umi_extracted/{sample}_umi_stats.txt")
            hist = read_family_histogram(
                f"{results_dir}/07_grouped/{sample}.family_size_histogram.txt"
            )
            report.append(
                f"{sample:<30} {umi['total_reads']:<18} "
                f"{umi['success_rate']:<12.1f} "
                f"{hist['total_families']:<14} "
                f"{hist['median_family_size']}"
            )

        report.append("")
        report.append("PIPELINE STEPS (PDF DOC-001337)")
        report.append("-" * 80)
        steps = [
            "1.  Custom  : UMI extracted from read-name header → prepended inline",
            "2.  Picard  : FastqToSam",
            "3.  fgbio   : ExtractUmisFromBam  (read-structure 5M+T, tags ZA/ZB/RX)",
            "4.  Picard  : SamToFastq (INTERLEAVE=true)",
            "5.  bwa mem : -p alignment",
            "6.  Picard  : MergeBamAlignment (CLIP_ADAPTERS=false, EXPECTED_ORIENTATIONS=FR)",
            "7.  fgbio   : GroupReadsByUmi --strategy=paired --raw-tag=RX --min-map-q=10 --edits=1",
            "8.  fgbio   : CallDuplexConsensusReads --error-rate-pre-umi=45 --min-reads 2 1 1",
            "9.  Picard  : SamToFastq (INTERLEAVE=true)",
            "10. bwa mem : -p re-alignment of consensus reads",
            "11. Picard  : MergeBamAlignment (retain ZS/ZI/ZM/ZC/ZN tags)",
        ]
        report.extend(steps)
        report.append("")

    with open(output_file, "w") as f:
        f.write("\n".join(report) + "\n")

    print("\n".join(report))


if __name__ == "__main__":
    import sys
    output_file = sys.argv[1] if len(sys.argv) > 1 else "results/summary_report.txt"
    generate_report(output_file)

This script collects statistics from all pipeline steps and generates
a human-readable summary report.
"""

import os
import json
import glob
from pathlib import Path


def read_umi_stats(stats_file):
    """Parse UMI extraction statistics."""
    stats = {
        'total_reads': 0,
        'umis_extracted': 0,
        'success_rate': 0
    }
    
    try:
        with open(stats_file, 'r') as f:
            for line in f:
                if 'Total reads processed:' in line:
                    stats['total_reads'] = int(line.split(':')[1].strip())
                elif 'UMIs successfully extracted:' in line:
                    stats['umis_extracted'] = int(line.split(':')[1].strip())
                elif 'Extraction success rate:' in line:
                    stats['success_rate'] = float(line.split(':')[1].strip().rstrip('%'))
    except Exception as e:
        print(f"Warning: Could not parse {stats_file}: {e}")
    
    return stats


def read_dedup_stats(stats_file):
    """Parse umi_tools deduplication statistics."""
    stats = {
        'reads_in': 0,
        'reads_out': 0,
        'duplicates_removed': 0
    }
    
    try:
        with open(stats_file, 'r') as f:
            content = f.read()
            # Parse JSON-like stats
            if 'reads_in' in content:
                import json
                data = json.loads(content)
                stats.update(data)
    except Exception as e:
        pass
    
    return stats


def read_trimming_stats(trimming_file):
    """Parse Trimmomatic trimming statistics."""
    stats = {
        'input_pairs': 0,
        'surviving_pairs': 0,
        'forward_only': 0,
        'reverse_only': 0
    }
    
    try:
        with open(trimming_file, 'r') as f:
            for line in f:
                if 'Input Read Pairs:' in line:
                    stats['input_pairs'] = int(line.split(':')[1].strip())
                elif 'Both Surviving:' in line:
                    stats['surviving_pairs'] = int(line.split(':')[1].strip())
                elif 'Forward Only Surviving:' in line:
                    stats['forward_only'] = int(line.split(':')[1].strip())
                elif 'Reverse Only Surviving:' in line:
                    stats['reverse_only'] = int(line.split(':')[1].strip())
    except Exception as e:
        pass
    
    return stats


def generate_report(output_file):
    """Generate comprehensive pipeline summary report."""
    
    # Determine results directory
    results_dir = "results"
    
    # Find all samples
    umi_stats_files = glob.glob(f"{results_dir}/02_umi_extracted/*_umi_stats.txt")
    samples = [Path(f).stem.replace('_umi_stats', '') for f in umi_stats_files]
    
    report = []
    report.append("=" * 80)
    report.append("TWIST PANEL UMI PREPROCESSING PIPELINE - SUMMARY REPORT")
    report.append("=" * 80)
    report.append("")
    
    if not samples:
        report.append("WARNING: No samples found in results directory")
        report.append("")
    else:
        # Sample processing summary
        report.append("SAMPLE PROCESSING SUMMARY")
        report.append("-" * 80)
        report.append(f"{'Sample':<30} {'Total Reads':<20} {'UMI Extr. %':<15} {'Final Reads':<15}")
        report.append("-" * 80)
        
        total_reads_all = 0
        total_umis_all = 0
        
        for sample in sorted(samples):
            # UMI stats
            umi_stats_file = f"{results_dir}/02_umi_extracted/{sample}_umi_stats.txt"
            umi_stats = read_umi_stats(umi_stats_file)
            
            # Trimming stats
            trimming_file = f"{results_dir}/03_trimmed/{sample}.trimming_report.txt"
            trim_stats = read_trimming_stats(trimming_file)
            
            # Dedup stats
            dedup_file = f"{results_dir}/05_dedup/{sample}_dedup_stats.txt"
            dedup_stats = read_dedup_stats(dedup_file)
            
            total_reads_all += umi_stats['total_reads']
            total_umis_all += umi_stats['umis_extracted']
            
            final_reads = trim_stats.get('surviving_pairs', 0)
            
            report.append(
                f"{sample:<30} {umi_stats['total_reads']:<20} "
                f"{umi_stats['success_rate']:<14.1f}% {final_reads:<15}"
            )
        
        report.append("-" * 80)
        report.append(f"{'TOTAL':<30} {total_reads_all:<20} "
                     f"{100 * total_umis_all / max(1, total_reads_all):<14.1f}%")
        report.append("")
        
        # Detailed stats by step
        report.append("")
        report.append("DETAILED PROCESSING STATISTICS BY STEP")
        report.append("=" * 80)
        report.append("")
        
        for sample in sorted(samples):
            report.append(f"Sample: {sample}")
            report.append("-" * 80)
            
            # Step 1: UMI Extraction
            umi_stats_file = f"{results_dir}/02_umi_extracted/{sample}_umi_stats.txt"
            umi_stats = read_umi_stats(umi_stats_file)
            report.append(f"  Step 1 - UMI Extraction from Headers:")
            report.append(f"    Input reads: {umi_stats['total_reads']}")
            report.append(f"    UMIs extracted: {umi_stats['umis_extracted']}")
            report.append(f"    Success rate: {umi_stats['success_rate']:.1f}%")
            report.append("")
            
            # Step 2: Trimming
            trimming_file = f"{results_dir}/03_trimmed/{sample}.trimming_report.txt"
            trim_stats = read_trimming_stats(trimming_file)
            report.append(f"  Step 2 - Adapter Trimming:")
            report.append(f"    Input pairs: {trim_stats['input_pairs']}")
            report.append(f"    Both surviving: {trim_stats['surviving_pairs']}")
            if trim_stats['input_pairs'] > 0:
                survival_rate = 100 * trim_stats['surviving_pairs'] / trim_stats['input_pairs']
                report.append(f"    Survival rate: {survival_rate:.1f}%")
            report.append("")
            
            # Step 3: Alignment (would need additional stats)
            report.append(f"  Step 3 - Reference Alignment:")
            bam_file = f"{results_dir}/04_aligned/{sample}.bam"
            if os.path.exists(bam_file):
                bam_size = os.path.getsize(bam_file) / (1024**3)  # Convert to GB
                report.append(f"    BAM file size: {bam_size:.2f} GB")
            report.append("")
            
            # Step 4: UMI Deduplication
            report.append(f"  Step 4 - UMI-Based Deduplication:")
            dedup_file = f"{results_dir}/05_dedup/{sample}_dedup_stats.txt"
            if os.path.exists(dedup_file):
                dedup_stats = read_dedup_stats(dedup_file)
                if dedup_stats['reads_in'] > 0:
                    report.append(f"    Input reads: {dedup_stats['reads_in']}")
                    report.append(f"    Output reads: {dedup_stats['reads_out']}")
                    dup_rate = 100 * (dedup_stats['reads_in'] - dedup_stats['reads_out']) / dedup_stats['reads_in']
                    report.append(f"    Duplicates removed: {dup_rate:.1f}%")
            else:
                report.append(f"    Deduplication statistics file not found")
            
            report.append("")
    
    # Pipeline parameters
    report.append("")
    report.append("PIPELINE CONFIGURATION")
    report.append("=" * 80)
    report.append("This pipeline follows the Twist Bioscience technical note:")
    report.append("  DOC-001337 - Processing Sequencing Data Utilizing UMI")
    report.append("")
    report.append("Key adaptations for pre-indexed UMIs:")
    report.append("  1. UMI extraction from fastq headers")
    report.append("  2. Prepending UMI sequences for downstream processing")
    report.append("  3. UMI-based deduplication via umi_tools")
    report.append("")
    
    # Output report
    with open(output_file, 'w') as f:
        f.write('\n'.join(report))
    
    # Also print to stdout
    print('\n'.join(report))


if __name__ == "__main__":
    import sys
    output_file = sys.argv[1] if len(sys.argv) > 1 else "results/summary_report.txt"
    generate_report(output_file)
