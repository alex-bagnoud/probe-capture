#!/usr/bin/env python3

import sys

blast_file = sys.argv[1]
kos_file = sys.argv[2]

ko_lookup = {}
with open(kos_file, 'r') as f:
    for line in f:
        line = line.strip()
        prot = line.split()[0]
        ko = line.split()[2]
        ko_lookup[prot] = ko

with open(blast_file, 'r') as f:
    for line in f:
        line = line.strip()
        frame = line.split()[0]
        prot = line.split()[1]
        try:
            ko = ko_lookup[prot]
        except KeyError:
            ko = "NA"

        print(f"{frame}\t{ko}")
