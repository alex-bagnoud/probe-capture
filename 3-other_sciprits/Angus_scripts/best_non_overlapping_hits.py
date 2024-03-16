#!/usr/bin/env python3

# gets up to two of the best hits if they are non-overlapping, per query. 
# otherwise just the single best hit is retrieved. Output file should follow
# the following dmnd parameterization:
# --outfmt 6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend\
#   sstart send evalue bitscore score positive gaps qlen slen
# Overlap tolerance defaults to 5 amino acids to account for possible error
# Output format should be sorted by quality of hits

import sys

OVLP_TOL = 5
dmnd_file = sys.argv[1]

def overlap(a, b):
    """a and b are tuples (or lists) of exactly 2 elements representing a
    range, returns the size of the overlap"""
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def get_range(line):
    return (int(line.split()[7]), int(line.split()[8]))


def insert(line, current_list):
    if overlap( get_range(line), get_range(current_list[0]) ) < OVLP_TOL:
        return [current_list[0], line]
    else:
        return current_list

def get_hits(fin):
    hits = {}
    for line in fin:
        line = line.strip()
        seq_id = line.split()[0]
        if seq_id not in hits:
            hits[seq_id] = [line]
        elif len(hits[seq_id]) < 2:
            hits[seq_id] = insert(line, hits[seq_id])

    for key in hits:
        for hit in hits[key]:
            print(hit)


with open(dmnd_file, 'r') as f:
    get_hits(f)

