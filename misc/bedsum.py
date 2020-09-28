#! /usr/bin/env python3

import sys

seqid2sum = {}

for line in sys.stdin:
    text = line.strip()
    if not text or text.startswith("#"):
        continue
    vals = text.split("\t")
    if vals[0] not in seqid2sum:
        seqid2sum[vals[0]] = 0
    # Assuming zero-based coords convention
    seqid2sum[vals[0]] += abs(int(vals[2]) - int(vals[1]))

for seqid in seqid2sum:
    print(f"{seqid}\t{seqid2sum[seqid]}")
