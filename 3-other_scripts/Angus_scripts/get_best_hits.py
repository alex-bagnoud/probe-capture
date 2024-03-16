#!/usr/bin/env python3
#
# gets the best hits from HMMs, including up to two non-overlapping hits.

import sys
h_file = sys.argv[1] 

def get_best(fin):
    best = {}
    for line in fin:
        if line[0] != '#':
            line = line.strip()
            fields = line.split()
            pr_name = fields[0]
            ko_name = fields[2]
            score = float(fields[5])
            e_value = float(fields[4])
            if pr_name not in best:
                best[pr_name] = line
            elif score > float(best[pr_name].split()[5]):
                best[pr_name] = line

    for key in best:
        line = best[key]
        fields = line.split()
        pr_name = fields[0]
        ko_name = fields[2]
        score = float(fields[5])
        e_value = float(fields[4])
        print(f"{pr_name}\t{ko_name}\t{score}\t{e_value}")

with open(h_file, 'r') as f:
    get_best(f)

