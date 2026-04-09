#!/usr/bin/env python3
"""Copy UMI values from read names into RX/ZA/ZB BAM tags.

This script expects UMIs to be present in the read name and uses a regex with
named groups 'umi_1' and 'umi_2' (same convention as extract_umi.py).

It streams BAM -> SAM -> BAM through samtools to avoid extra Python deps.
"""

import argparse
import re
import subprocess
import sys


def parse_args():
    parser = argparse.ArgumentParser(description="Annotate BAM with UMI tags from QNAME")
    parser.add_argument("--input-bam", required=True, help="Input BAM path")
    parser.add_argument("--output-bam", required=True, help="Output BAM path")
    parser.add_argument(
        "--umi-pattern",
        default=r"(?P<umi_1>[ACGTN]+)\+(?P<umi_2>[ACGTN]+)",
        help="Regex with named groups 'umi_1' and 'umi_2'",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    pattern = re.compile(args.umi_pattern)

    in_cmd = ["samtools", "view", "-h", args.input_bam]
    out_cmd = ["samtools", "view", "-b", "-o", args.output_bam, "-"]

    p_in = subprocess.Popen(in_cmd, stdout=subprocess.PIPE, text=True)
    p_out = subprocess.Popen(out_cmd, stdin=subprocess.PIPE, text=True)

    total = 0
    tagged = 0
    failed = 0

    try:
        assert p_in.stdout is not None
        assert p_out.stdin is not None

        for line in p_in.stdout:
            if line.startswith("@"):
                p_out.stdin.write(line)
                continue

            total += 1
            fields = line.rstrip("\n").split("\t")
            qname = fields[0]
            match = pattern.search(qname)

            if not match:
                failed += 1
                if failed <= 5:
                    print(f"No UMI match for read name: {qname}", file=sys.stderr)
                p_out.stdin.write(line)
                continue

            gd = match.groupdict()
            umi_1 = gd.get("umi_1")
            umi_2 = gd.get("umi_2")
            if not umi_1 or not umi_2:
                failed += 1
                if failed <= 5:
                    print(f"UMI regex missing named groups for read: {qname}", file=sys.stderr)
                p_out.stdin.write(line)
                continue

            tags = [t for t in fields[11:] if not (t.startswith("RX:Z:") or t.startswith("ZA:Z:") or t.startswith("ZB:Z:"))]
            tags.extend([
                f"ZA:Z:{umi_1}",
                f"ZB:Z:{umi_2}",
                f"RX:Z:{umi_1}-{umi_2}",
            ])

            p_out.stdin.write("\t".join(fields[:11] + tags) + "\n")
            tagged += 1

    finally:
        if p_in.stdout is not None:
            p_in.stdout.close()
        if p_out.stdin is not None:
            p_out.stdin.close()

    rc_in = p_in.wait()
    rc_out = p_out.wait()
    if rc_in != 0 or rc_out != 0:
        raise SystemExit(f"samtools pipeline failed (in={rc_in}, out={rc_out})")

    print(f"Processed reads: {total}")
    print(f"Tagged reads:    {tagged}")
    print(f"Unmatched reads: {failed}")

    if failed > 0:
        raise SystemExit(
            "UMI tagging failed for one or more reads. Check --umi-pattern and read name format."
        )


if __name__ == "__main__":
    main()
