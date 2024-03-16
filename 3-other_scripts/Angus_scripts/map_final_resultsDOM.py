#!/usr/bin/env python3 

import sys
read_to_pr_fname = sys.argv[1]  # i.e., top2 blast results
pr_to_ko_fname = sys.argv[2]    # i.e., hmm results for genomes
read_to_ko_fname = sys.argv[3]  # i.e., actual results/best HMM per frame

def build_read_to_pr(fin):
    r_to_p = {}
    for line in fin:
        line = line.strip()
        read = line.split()[0]
        prot = line.split()[1]
        r_to_p[read] = prot
    return r_to_p

def build_pr_to_ko(fin):
    p_to_k = {}
    for line in fin:
        line = line.strip()
        pr = line.split()[0]
        ko = line.split()[2]
        p_to_k[pr] = ko
    return p_to_k

def get_results(fin, p_to_k, r_to_p):
    print(f"#read\tpr\te_val\tbitscore\tdom_score\tcovg\tpred_ko\ttrue_ko")
    for line in fin:
        line = line.strip()
        fields = line.split()
        read = fields[0]
        e_val = float(fields[6])
        bitscore = float(fields[7])
        dom_score = float(fields[13])
        pred_ko = fields[3]
        length = int(fields[2])
        aln_from = int(fields[17])
        aln_to = int(fields[18])
        covg = float(aln_to - aln_from + 1)/length
        true_ko = "NA"
        pr = "NA"
        try:
            pr = r_to_p[read]
        except KeyError:
            pr = "NA"
        try:
            true_ko = p_to_k[pr]
        except KeyError:
            true_ko = "NA"

        print(f"{read}\t{pr}\t{e_val}\t{bitscore}\t{dom_score}\t{covg}\t{pred_ko}\t{true_ko}")

read_to_pr = {}
pr_to_ko = {}
with open(read_to_pr_fname, 'r') as f:
    read_to_pr = build_read_to_pr(f)
with open(pr_to_ko_fname, 'r') as f:
    pr_to_ko = build_pr_to_ko(f)
with open(read_to_ko_fname, 'r') as f:
    get_results(f, pr_to_ko, read_to_pr)
    
