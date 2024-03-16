#!/usr/bin/env python3

# input list of files in format gene_...extension
# generates stats from the input files
# Note, second column should be the correct KO
# File 1 should have the true positive mappings (read followed by KO list)
# File 2 should have all positive counts

import sys
ko_file = sys.argv[1]
counts_file = sys.argv[2]

all_pos = {}
with open(counts_file, 'r') as f:
    for line in f:
        line = line.strip()
        ko = line.split()[1]
        count = int(line.split()[2])
        all_pos[ko] = count

read_mapping = {}
with open(ko_file, 'r') as f:
    for line in f:
        line = line.strip()
        read = line.split()[0]
        kos = line.split()[1].split(",")
        read_mapping[read] = kos

for f in sys.argv[3:]:
    with open(f, 'r') as fin:
        tp_count = 0
        fp_count = 0
        pred_ko = ""
        for line in fin:
            line = line.strip()
            read = line.split()[0]
            pred_ko = line.split()[1]
            true_ko = line.split()[2]
            try:
                if pred_ko == true_ko:
                    tp_count = tp_count + 1
                else:
                    fp_count = fp_count + 1
            except KeyError:
                fp_count = fp_count + 1
        true_count = all_pos[pred_ko]
        fn_count = true_count - tp_count
        precision = float(tp_count)/(tp_count + fp_count)
        recall = float(tp_count)/(tp_count + fn_count)
        gene = fin.name.split("_")[0].strip()
        print(f"{gene}\t{tp_count}\t{fp_count}\t{fn_count}\t{precision:.4f}\t{recall:.4f}")





