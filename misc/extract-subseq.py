#! /usr/bin/env python3

import sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument(
    "--regions", default=False, help="bedfile of regions to extract, 0-based indices"
)
parser.add_argument(
    "--one-based-coords",
    action="store_true",
    default=False,
    help="the REGIONS bedfile uses 1-based coordinates rather than 0-based",
)
args = parser.parse_args()

if not args.regions:
    exit(0)

if args.one_based_coords:
    offset = 1
else:
    offset = 0

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
        score = region[3]
        strand = region[4].strip()
        if not strand:
            strand = "+"

        if strand == "+" or strand == ".":
            seqtext = str(r.seq[start - offset : end])
            print(
                f'>{locus} source="{chr}" start{offset}="{start}" end{offset}="{end}" strand="{strand}"'
            )
            print(f"{seqtext}")

        elif strand == "-":
            seqtext = str(r.seq[start - offset : end].reverse_complement())
            print(
                f'>{locus} source="{chr}" start{offset}="{start}" end{offset}="{end}" strand="{strand}"'
            )
            print(f"{seqtext}")
            ##print(f"{str(r.seq[start : end])}")
