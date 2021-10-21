#! /usr/bin/env python

import argparse
import sys
from Bio import SeqIO


def get_options():
    p = argparse.ArgumentParser()
    p.add_argument(
        "infile", nargs="?", default=sys.stdin, help="Filename, default stdin."
    )
    p.add_argument(
        "-l",
        "--len",
        type=int,
        default=0,
        help="Produce reads this long when possible, default (0) full length.",
    )
    return p.parse_args()


if __name__ == "__main__":
    args = get_options()

    targ_len = args.len

    for r in SeqIO.parse(args.infile, "fastq"):
        seq_len = len(r.seq)
        start = 0
        end = seq_len
        if targ_len > 0 and targ_len < seq_len:
            d = int((seq_len - targ_len) / 2)
            start = d
            end = start + targ_len
        ##print(f"Sequence length is {seq_len}. The fillet is {start} to {end}.")
        try:
            print(r[start:end].format("fastq"))
        except:
            print(f"fq-fillet.bash failed on sequence {r.id}, exiting") 
            exit(1)
