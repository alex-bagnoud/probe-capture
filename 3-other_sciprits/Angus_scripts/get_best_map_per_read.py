#!/usr/bin/env python3

import sys
with open(sys.argv[1], 'r') as f:
    reads = {}
    for line in f:
        line = line.strip()
        read = line.split()[0]#[:-2]
        if read not in reads:
            reads[read] = line
        elif float(line.split()[11]) < float(reads[read].split()[11]):
            reads[read] = line

    for key in reads:
        print(reads[key])
