#! /usr/bin/env python3

import numpy as np

ngenes = 30
nsamps_a = 2
nsamps_b = 3
header = []
header.append("gene")
header.extend([f"sample{1 + int(s) :02}" for s in range(nsamps_a + nsamps_b)])
print("\t".join(header))
for i in range(ngenes):
    record = []
    record.append(f"gene{1 + i :03}")
    record.extend([str(x) for x in np.random.poisson(lam=100, size=nsamps_a)])
    record.extend([str(x) for x in np.random.poisson(lam=100, size=nsamps_b)])
    print("\t".join(record))
