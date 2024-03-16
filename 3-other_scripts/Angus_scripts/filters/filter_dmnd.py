#!/usr/bin/env python3
# Filter diamond results based on 95.0% identity and print only the best hit for each query sequence

import sys
table = sys.argv[1]

def filt(fin):
    results = {}

    for line in fin:
        line = line.strip()
        fields = line.split()
        qid = fields[0]
        sid = fields[1]
        ident = float(fields[2])

        if ident > 95.0:
            if qid not in results:
                results[qid] = (sid, ident)
            elif ident > results[qid][1]:
                results[qid] = (sid, ident)

    for key in results:
        print(f"{key}\t{results[key][0]}\t{results[key][1]}")


with open(table, 'r') as f:
    filt(f)
