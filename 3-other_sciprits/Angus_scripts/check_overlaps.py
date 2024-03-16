#!/usr/bin/env python3

# checks for potential protein fusions

import sys
dom_file = sys.argv[1]

def overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def insert(line, current):
    e_val = float(line.split()[6])

    if len(current == 0):
        return [line]
    elif len(current == 1):
        if e_val < float(current[0].split()[6]):
            return [line, current[0]]
        else:
            return [current[0], line]
    
    if e_val < float(current[0].split()[6]):
        return [line, current[0]]
    elif e_val < float(current[1].split()[6]):
        return [current[0], line] 

def check_overlaps(fin):
    hits = {}
    for line in fin:
        line = line.strip()
        prot = line.split()[0]
        if prot not in hits:
            hits[prot] = [line]
        else:
            hits[prot] = insert(line, hits[prot])

    for key in hits:
        if len(hits[key] == 2):
            c1 = int(hits[key][0].split()[17])
            c2 = int(hits[key][0].split()[18])
            c3 = int(hits[key][1].split()[17])
            c4 = int(hits[key][1].split()[18])
            if overlap( (c1,c2), (c3,c4) ) <= 10:
                print(hits[key])
