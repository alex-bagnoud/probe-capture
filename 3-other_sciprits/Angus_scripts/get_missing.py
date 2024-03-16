#!/usr/bin/env python3
import sys

def get_missing(fin):
    missing = set()
    for line in fin:
        if int(line.split()[1]) == 0:
            missing.add(line.split()[0].strip())
    return missing

def print_missing(fin, missing):
    for line in fin:
        ident = line.split()[0][:-2]
        if ident in missing:
            print(line.strip())

miss_set = set()
with open(sys.argv[1], 'r') as f:
    miss_set = get_missing(f)
with open(sys.argv[2], 'r') as f:
    print_missing(f, miss_set)
    
