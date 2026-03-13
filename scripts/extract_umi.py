#!/usr/bin/env python3
"""
Extract UMI sequences from fastq read-name headers and prepare reads for
the fgbio FastqToBam workflow.

Input fastq header format:
  @VH01088:45:AAJ2V5CM5:1:1101:60230:1000:GTCGT+CAGGA 1:N:0:AATGAGGTAA+GTATAACCTG
                                            ^^^^^  ^^^^^
                                            UMI1   UMI2  (5 bp each by default)

Two modes:
  Default (no --split_umis): full "UMI1+UMI2" string prepended to R1 only.
  --split_umis (required for fgbio): UMI1 prepended to R1, UMI2 prepended to
    R2.  fgbio FastqToBam then strips them via --read-structures <N>M+T <N>M+T.
"""

import argparse
import gzip
import re
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract UMIs from fastq headers and prepend inline"
    )
    parser.add_argument("--r1", required=True, help="Input R1 fastq (.fastq/.fastq.gz)")
    parser.add_argument("--r2", required=True, help="Input R2 fastq (.fastq/.fastq.gz)")
    parser.add_argument("--out_r1", required=True, help="Output R1 fastq")
    parser.add_argument("--out_r2", required=True, help="Output R2 fastq")
    parser.add_argument("--stats", required=True, help="Output statistics file")
    parser.add_argument(
        "--umi_pattern",
        default=r"(?P<umi_1>[ACGTN]+)\+(?P<umi_2>[ACGTN]+)",
        help=(
            "Regex with named groups 'umi_1' and 'umi_2'.  "
            "Default captures any two sequences separated by '+'."
        ),
    )
    parser.add_argument(
        "--split_umis",
        action="store_true",
        default=False,
        help=(
            "Split UMI1/UMI2 and prepend each to its own read "
            "(required for fgbio FastqToBam --read-structures).  "
            "Without this flag the full 'UMI1+UMI2' string is prepended to R1 only."
        ),
    )
    return parser.parse_args()


def open_fastq(filename):
    return gzip.open(filename, "rt") if filename.endswith(".gz") else open(filename)


def open_output(filename):
    return gzip.open(filename, "wt") if filename.endswith(".gz") else open(filename, "w")


def extract_umis(header, compiled_pattern):
    """
    Return (umi_1, umi_2) strings from the read name, or (None, None).
    Named groups 'umi_1' and 'umi_2' in the regex are required.
    """
    m = compiled_pattern.search(header)
    if not m:
        return None, None
    gd = m.groupdict()
    return gd.get("umi_1"), gd.get("umi_2")


def annotate_header(header, umi_1, umi_2):
    """
    Append ':umi:<UMI1>+<UMI2>' to the instrument field of the read name.
    The comment (part after the space) is left unchanged.
    """
    header = header.rstrip()
    parts = header.split(" ", 1)
    parts[0] = f"{parts[0]}:umi:{umi_1}+{umi_2}"
    return " ".join(parts)


def process_fastq_pair(r1_file, r2_file, out_r1, out_r2, compiled_pattern, split_umis):
    """
    Process paired-end fastq files.

    split_umis=True  → UMI1 prepended to R1, UMI2 prepended to R2
                        (needed for fgbio --read-structures 5M+T 5M+T)
    split_umis=False → full UMI1+UMI2 prepended to R1 only
    """
    stats = {
        "total_reads": 0,
        "umis_extracted": 0,
        "umis_failed": 0,
        "umi1_sequences": defaultdict(int),
        "umi2_sequences": defaultdict(int),
        "errors": [],
    }

    with (
        open_fastq(r1_file) as f_r1,
        open_fastq(r2_file) as f_r2,
        open_output(out_r1) as o_r1,
        open_output(out_r2) as o_r2,
    ):
        while True:
            h1 = f_r1.readline()
            if not h1:
                break
            s1, p1, q1 = f_r1.readline(), f_r1.readline(), f_r1.readline()
            h2, s2, p2, q2 = (
                f_r2.readline(), f_r2.readline(),
                f_r2.readline(), f_r2.readline(),
            )

            stats["total_reads"] += 1
            umi_1, umi_2 = extract_umis(h1, compiled_pattern)

            if umi_1 and umi_2:
                stats["umis_extracted"] += 1
                stats["umi1_sequences"][umi_1] += 1
                stats["umi2_sequences"][umi_2] += 1

                new_h1 = annotate_header(h1, umi_1, umi_2) + "\n"
                new_h2 = annotate_header(h2, umi_1, umi_2) + "\n"

                if split_umis:
                    # UMI1 → R1 prefix, UMI2 → R2 prefix
                    # Quality for the added UMI bases is set to max (I = Phred 40)
                    new_s1 = umi_1 + s1.rstrip() + "\n"
                    new_q1 = ("I" * len(umi_1)) + q1.rstrip() + "\n"
                    new_s2 = umi_2 + s2.rstrip() + "\n"
                    new_q2 = ("I" * len(umi_2)) + q2.rstrip() + "\n"
                else:
                    combined = umi_1 + umi_2
                    new_s1 = combined + s1.rstrip() + "\n"
                    new_q1 = ("I" * len(combined)) + q1.rstrip() + "\n"
                    new_s2, new_q2 = s2, q2

                o_r1.write(new_h1); o_r1.write(new_s1); o_r1.write(p1); o_r1.write(new_q1)
                o_r2.write(new_h2); o_r2.write(new_s2); o_r2.write(p2); o_r2.write(new_q2)

            else:
                stats["umis_failed"] += 1
                if len(stats["errors"]) < 20:
                    stats["errors"].append(h1.rstrip())
                o_r1.write(h1); o_r1.write(s1); o_r1.write(p1); o_r1.write(q1)
                o_r2.write(h2); o_r2.write(s2); o_r2.write(p2); o_r2.write(q2)

            if stats["total_reads"] % 100_000 == 0:
                print(f"  Processed {stats['total_reads']:,} read pairs ...", flush=True)

    return stats


def write_statistics(stats, output_file):
    total = stats["total_reads"]
    with open(output_file, "w") as f:
        f.write("UMI Extraction Statistics\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total reads processed:      {total}\n")
        f.write(f"UMIs successfully extracted: {stats['umis_extracted']}\n")
        f.write(f"UMIs failed to extract:      {stats['umis_failed']}\n")
        f.write(
            f"Extraction success rate:     "
            f"{100 * stats['umis_extracted'] / max(1, total):.2f}%\n\n"
        )

        for label, key in [("UMI1", "umi1_sequences"), ("UMI2", "umi2_sequences")]:
            f.write(f"Top 20 {label} sequences:\n")
            for seq, cnt in sorted(stats[key].items(), key=lambda x: -x[1])[:20]:
                f.write(f"  {seq}: {cnt}\n")
            f.write("\n")

        if stats["errors"]:
            f.write("Headers that failed UMI extraction (first 20):\n")
            for h in stats["errors"]:
                f.write(f"  {h}\n")


if __name__ == "__main__":
    args = parse_args()

    print("Extracting UMIs from fastq headers ...")
    print(f"  Input R1      : {args.r1}")
    print(f"  Input R2      : {args.r2}")
    print(f"  UMI pattern   : {args.umi_pattern}")
    print(f"  Split UMIs    : {args.split_umis}")
    print()

    compiled = re.compile(args.umi_pattern)
    stats = process_fastq_pair(
        args.r1, args.r2,
        args.out_r1, args.out_r2,
        compiled,
        split_umis=args.split_umis,
    )

    write_statistics(stats, args.stats)

    total = stats["total_reads"]
    print(f"\nDone.  {stats['umis_extracted']:,}/{total:,} reads had UMIs extracted "
          f"({100 * stats['umis_extracted'] / max(1, total):.1f}%)")
