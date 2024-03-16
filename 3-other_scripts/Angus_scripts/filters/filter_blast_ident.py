#!/usr/bin/env python3
# filter dmnd blast results (argv 1) based on identity (argv 2)

import sys

bname = sys.argv[1]
filt = float(sys.argv[2])

def filter(fin, f_value=filt):
    for line in fin:
        ident = float(line.split()[2])
        if ident >= f_value:
            print(line.strip())


with open(bname, 'r') as f:
    filter(f)
