#!/usr/bin/env python3

import sys

fname = sys.argv[1]
ko = sys.argv[2]

with open(sys.argv[1], 'r') as fin:
    for line in fin:
        line = line.strip()
        print(f"{line}\t{ko}")
