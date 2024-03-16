#!/usr/bin/env python3

# get the top KO mapping (by e-value) and print for each protein

import sys

with open(sys.argv[1], 'r') as fin:
    mapping = {}
    for line in fin:
        line = line.strip()
        if line[0] != "#":
            pr_id = line.split()[0]
            e_val = float(line.split()[4])
            if pr_id not in mapping:
                mapping[pr_id] = line
            elif e_val < float(mapping[pr_id].split()[4]):
                mapping[pr_id] = line
    for key in mapping:
        print(mapping[key])
