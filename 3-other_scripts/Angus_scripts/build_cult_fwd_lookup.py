#!/usr/bin/env python3

import sys

blast_file = sys.argv[1]
kofam_file = sys.argv[2]

def build_kofam_lookup(fin):
    lookup = {}
    for line in fin:
        line = line.strip()
        prot = line.split()[0]
        ko = line.split()[1]
        lookup[prot] = ko
    return lookup

def get_best_prot(fin):
    read_to_pr = {}
    for line in fin:
        line = line.strip()
        if int(line.split()[0][-1]) not in [1,2,3]:
            continue
        seq_id = line.split()[0][:-2]
        prot_id = line.split()[1]
        evalue = float(line.split()[11])
        if seq_id not in read_to_pr:
            read_to_pr[seq_id] = (evalue, prot_id)
        elif evalue < read_to_pr[seq_id][0]:
            read_to_pr[seq_id] = (evalue, prot_id)

    return read_to_pr

def print_res(rtp, pup):
    for key in rtp:
        prot = rtp[key][1]
        if prot in pup:
            ko = pup[prot]
        else:
            ko = "NA"
        print(f"{key}\t{ko}")

p = {}
r = {}
with open(kofam_file, 'r') as f:
    p = build_kofam_lookup(f)
with open(blast_file, 'r') as f:
    r = get_best_prot(f)
print_res(r, p)
