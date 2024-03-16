#!/usr/bin/env python3

# short script to remove identical reads using dictionaries
# Default is to output FASTA format, output table with mappings
# specified with second arg.

# NOTE: The 2 line FASTA file assumption is reasonable here.
# But anyway, this script will verify that

import sys

fastaname = sys.argv[1]
try:
    outtable = sys.argv[2]
except IndexError:
    outtable = None

def dedup_results(fin):
    results = {}
    while True:
        header = fin.readline().strip()
        seq = fin.readline().strip()
    
        if (not header) or (not seq): break

        if header[0] != ">":
            raise TypeError

        if seq not in results:
            results[seq] = []
        results[seq].append(header)
    return results


def print_results(deduped, output=None):
    for key in deduped:
        print(deduped[key][0])
        print(key)

        if output:
            output.write("\t".join(deduped[key]))

res = {}
with open(fastaname, 'r') as fhandle:
    res = dedup_results(fhandle)

if outtable:
    with open(outtable, 'w') as fhandle:
        print_results(res, fhandle)
else:
    print_results(res)
