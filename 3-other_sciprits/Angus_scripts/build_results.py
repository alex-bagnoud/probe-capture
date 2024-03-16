#!/usr/bin/env python3
import sys

genome_ko = sys.argv[1]     # ko mappings for the genome files
results_ko = sys.argv[2]    # ko hmm file for samples
mapfile = sys.argv[3]       # read mapping for samples

def build_prot_lookup(fin):
    prot_lookup = {}
    for line in fin:
        line = line.strip()
        fields = line.split()
        try:
            prot_lookup[fields[0]] = fields[1:]
        except IndexError:
            prot_lookup[fields[0]] = ["NA"]*4
    return prot_lookup

def build_mapping(fin):
    mapping = {}
    for line in fin:
        line = line.strip()
        fields = line.split()
        read = fields[0]
        prot = fields[1]
        mapping[read] = prot
    return mapping

def build_results(fin, prot_ko, mapping):
    for line in fin:
        if line[0] != '#':
            if line[0] == "*":
                line = line[1:]

            fields = line.strip().split()
            read = fields[0]
            KO = fields[1]
            bitscore = fields[2]
            e_val = fields[3]

            pdata = "\t".join(["NA"]*3)
            if read in mapping and mapping[read] in prot_ko:
                pdata = "\t".join(prot_ko[mapping[read]])

            print(f"{read}\t{bitscore}\t{e_val}\t{KO}\t{pdata}")

p_lookup = {}
mapping_lookup = {}

with open(genome_ko, 'r') as fin:
    p_lookup = build_prot_lookup(fin)
with open(mapfile, 'r') as fin:
    mapping_lookup = build_mapping(fin)
with open(results_ko, 'r') as fin:
    build_results(fin, p_lookup, mapping_lookup)
