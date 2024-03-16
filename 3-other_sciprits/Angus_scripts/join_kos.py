#!/usr/bin/env python3

import sys
with open(sys.argv[1], 'r') as f:
    ko_lookup = {}
    for line in f:
        line = line.strip()
        read = line.split()[0][:-2]
        ko = line.split()[1]
        if read not in ko_lookup:
            ko_lookup[read] = set()
        ko_lookup[read].add(ko)
    
    for key in ko_lookup:
        kos = ",".join(list(ko_lookup[key]))
        print(f"{key}\t{kos}")

