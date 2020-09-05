#! /usr/bin/env python3

import sys
import argparse
from Bio import SeqIO

# 23456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789
parser = argparse.ArgumentParser(
    description="Extract subsequences associated with specified regions in a "
    "reference.",
    epilog="A set of reference sequences is provided in FASTA format on stdin "
    "and regions of interest are specified in a bedfile. The subsequence "
    "associated with each region is written to stdout. When a strand is specified "
    "for a region (col 6 of the bedfile) the subsequence is reported on that "
    "strand.",
)
parser.add_argument(
    "-r",
    "--regions",
    default=False,
    help="bedfile of regions to extract, 0-based indices",
)
parser.add_argument(
    "-1",
    "--one-based-coords",
    action="store_const",
    const=1,
    dest="offset",
    default=0,
    help="the REGIONS bedfile uses 1-based coordinates rather than standard 0-based",
)
args = parser.parse_args()

if not args.regions:
    exit(0)

regions_by_chr = {}

with open(args.regions, "r") as infh:
    for line in infh:
        text = line.strip()
        if not text or text.isspace() or text.startswith("#"):
            continue
        vals = text.split("\t")
        ##print(vals)
        if vals[0] not in regions_by_chr:
            regions_by_chr[vals[0]] = []
        regions_by_chr[vals[0]].append(vals[1:])

for r in SeqIO.parse(sys.stdin, "fasta"):
    chr = r.id

    if chr not in regions_by_chr:
        continue

    for region in regions_by_chr[chr]:
        start = int(region[0])
        end = int(region[1])
        locus = region[2].strip()
        score = region[3]  # ignored
        strand = region[4].strip()

        if not strand or strand == "+" or strand == ".":
            seqtext = str(r.seq[start - args.offset : end])
        elif strand == "-":
            seqtext = str(r.seq[start - args.offset : end].reverse_complement())
        else:
            print(
                f"[WARNING] Skipping region {chr}\t{start}\t{end}\t{locus}\t{score}\t{strand}",
                file=sys.stderr,
            )
            continue
        print(
            f'>{locus} source="{chr}" start{args.offset}="{start}" end{args.offset}="{end}" strand="{strand}"'
        )
        print(f"{seqtext}")
