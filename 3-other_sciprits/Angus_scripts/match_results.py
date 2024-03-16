#!/usr/bin/env python3

import sys
mapping_file = sys.argv[1] # reads to correct ko
predicted_file = sys.argv[2] # HA method

def build_lookup(fin):
    lookup = {}
    for line in fin:
        line = line.strip()
        seq_id = line.split()[0]
        ko = line.split()[1]
        lookup[seq_id] = ko
    return lookup

def match_res(fin, lookup):
    for line in fin:
        line = line.strip()
        seq_id = line.split()[0]
        pred_ko = line.split()[1]
        try:
            true_ko = lookup[seq_id]
        except KeyError:
            true_ko = "NA"
        print(f"{seq_id}\t{pred_ko}\t{true_ko}")

lu = {}
with open(mapping_file, 'r') as f:
    lu = build_lookup(f)
with open(predicted_file, 'r') as f:
    match_res(f, lu)
