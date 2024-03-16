#!/usr/bin/env python3 

import sys

bname = sys.argv[1]
qname = sys.argv[2]

def build_qdict(fin):
    qd = {}
    for line in fin:
        if line[0] == ">":
            ident = line.strip().split()[0][1:-2]
            qd[ident] = 0
    return qd


def count_res(fin, qd):
    for line in fin:
        ident = line.strip().split()[0][:-2]
        qd[ident] = qd[ident] + 1

    for key in qd:
        print(f"{key}\t{qd[key]}")

qdict = {}
with open(qname, 'r') as f:
    qdict = build_qdict(f)
with open(bname, 'r') as f:
    count_res(f, qdict)
