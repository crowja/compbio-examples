#! /usr/bin/env python3

import argparse
import random
import sys

def get_command_line():
    ap = argparse.ArgumentParser()
    ap.add_argument("--num-indeps", type=int, default=1, help="Number of independent variables.")
    ap.add_argument("--num-deps", type=int, default=1, help="Number of dependent variables.")
    ap.add_argument("--num-recs", type=int, default=1, help="Number of records.")
    return ap.parse_args()

if __name__ == "__main__":
    args = get_command_line()

    # Print the header
    header = []
    header.append("LABEL")
    for j in range(args.num_indeps):
        header.append(f"INDEP_{j + 1 :03}")
    for j in range(args.num_deps):
        header.append(f"DEP_{j + 1 :03}")
    print(",".join(header))

    # Print the records
    for i in range(args.num_recs):
        record = []
        x = random.choice(['control', 'treated'])
        record.append(x)
        for j in range(args.num_indeps):
            x = random.choice(['a', 'b', 'c'])
            record.append(f"{j + 1}_{x}")
        for j in range(args.num_deps):
            x = random.random() * 100
            record.append(f"{x :0.4}")
        print(",".join(record))

