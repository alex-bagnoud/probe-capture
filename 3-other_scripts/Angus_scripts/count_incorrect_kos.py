#!/usr/bin/env python3

import sys
ko_file = sys.argv[1]
fname = sys.argv[2]
# checks if last col KO matches first, if not, adds last call KO to counts
# for all kos

ko_names = {}
with open(ko_file, 'r') as fin:
    for line in fin:
        ko, name = tuple(line.strip().split("\t"))
        ko = ko.split(":")[1]
        ko_names[ko] = name

counts = {}
with open(fname, 'r') as fin:
    counter = 0
    for line in fin:
        line = line.strip()
        ko1 = line.split()[-2]
        ko2 = line.split()[-1]
        if ko1 != ko2:
            counter = counter + 1
            if ko2 not in counts:
                counts[ko2] = 0
            counts[ko2] = counts[ko2] + 1

    for key in counts:
        prop = 100 * counts[key]/float(counter)
        if prop > 0.5:
            try:
                ko_name = ko_names[key]
            except KeyError:
                ko_name = ""
            print(f"{key}\t{counts[key]}\t{prop:.1f}\t{ko_name}")




