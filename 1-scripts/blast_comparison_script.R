setwd("~/OneDrive/Work/Collaboration/henri_probe_capture/5-blast_comparison/")

```{bash}

# Prepare QC-filtered amplicons reads for the BLAST analysis
cat ../1-amplicon_seq/temp2-qc_filtered_fasta_files/H2-1*_R1.fasta > 1a-qc_r1.fasta
cat ../1-amplicon_seq/temp2-qc_filtered_fasta_files/H2-1*_R2.fasta > 1b-qc_r2.fasta

usearch_v8.1.1861 --derep_fulllength 1a-qc_r1.fasta --fastaout 2a-r1_derep.fasta --sizeout
# 91091 seqs, 22975 uniques, 17905 singletons (77.9%)

usearch_v8.1.1861 --derep_fulllength 1b-qc_r2.fasta --fastaout 2b-r2_derep.fasta --sizeout
# 91905 seqs, 22004 uniques, 17213 singletons (78.2%)

# Create BLAST databases
makeblastdb -in 2a-r1_derep.fasta -dbtype nucl -out 3a-r1
makeblastdb -in 2b-r2_derep.fasta -dbtype nucl -out 3b-r2
makeblastdb -in ../2-probe_capture/1-derep_pc.fasta -dbtype nucl -out 3c-pc

# Run BLAST analysis
blastn -query 2a-r1_derep.fasta -db 3c-pc -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send slen evalue bitscore" \
-out 4a-r1_vs_pc_tab.txt
blastn -query 2a-r1_derep.fasta -db 3c-pc -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send slen evalue bitscore" \
-out 4a-r1_vs_pc_max5.txt -max_target_seqs 5
blastn -query 2a-r1_derep.fasta -db 3c-pc -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send slen evalue bitscore" \
-out 4a-r1_vs_pc_max1.txt -max_target_seqs 1


blastn -query 2b-r2_derep.fasta -db 3c-pc -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send slen evalue bitscore" \
-out 4b-r2_vs_pc_tab_max1.txt -max_target_seqs 1

blastn -query ../2-probe_capture/1-derep_pc.fasta -db 3a-r1 -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send slen evalue bitscore" \
-out 5a-pc_vs_r1_tab_max1.txt -max_target_seqs 1

blastn -query ../2-probe_capture/1-derep_pc.fasta -db 3a-r2 -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send slen evalue bitscore" \
-out 5b-pc_vs_r2_tab_max1.txt -max_target_seqs 1
```

r1.pc.blast <- read.table("4a-r1_vs_pc_tab.txt", sep = '\t')

# Keep the best hit for each query
r1.pc.blast.1seq <- r1.pc.blast[!duplicated(r1.pc.blast$V1),]
rm(r1.pc.blast)

# Name columns
colnames(r1.pc.blast.1seq) <- c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "slen", "evalue", "bitscore")

hist(r1.pc.blast.1seq$pident, breaks = seq(86, 100, by = 1))

# Compute 'final scores', which are alignement length * pident / query length
r1.pc.blast.1seq$fscore <- r1.pc.blast.1seq$pident * r1.pc.blast.1seq$length / 200
hist(r1.pc.blast.1seq$fscore, breaks = seq(44, 101, by = 1))


# Select alignement with pident >= 99 and no gap open
r1.pc.blast.1seq.99pident <- r1.pc.blast.1seq[r1.pc.blast.1seq$pident >= 99 & r1.pc.blast.1seq$gapopen == 0,]
dim(r1.pc.blast.1seq.99pident)
# 17805    14
table(r1.pc.blast.1seq.99pident$fscore)

# How many of those have a f-score bigger than 99% ?
sum(r1.pc.blast.1seq.99pident$fscore >= 99)
# 17012

# Total number of original sequences
nrow(r1.pc.blast.1seq)
# 22975

# Identify sequence with a lower fscore because uncomplete alignement
others <- r1.pc.blast.1seq.99pident[r1.pc.blast.1seq.99pident$fscore < 99,]

sum(others$send == 1)
# 71

sum (others$sstart == others$slen)
# 35

## Compare these results with max target seq = 1

r1.pc.blast.mts1 <- read.table("4a-r1_vs_pc_max1.txt", sep = '\t')

# Name columns
colnames(r1.pc.blast.mts1 ) <- c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "slen", "evalue", "bitscore")

hist(r1.pc.blast.mts1$pident, breaks = seq(86, 100, by = 1))
sum(r1.pc.blast.mts1$pident >= 99)

# Compute 'final scores', which are alignement length * pident / query length
r1.pc.blast.mts1$fscore <- r1.pc.blast.mts1$pident * r1.pc.blast.mts1$length / 200
hist(r1.pc.blast.mts1$fscore, breaks = seq(44, 101, by = 1))
sum(r1.pc.blast.mts1$fscore >= 99)

# Select alignement with pident >= 99 and no gap open
r1.pc.blast.mts1.99pident <- r1.pc.blast.mts1[r1.pc.blast.mts1$pident >= 99 & r1.pc.blast.mts1$gapopen == 0,]
dim(r1.pc.blast.mts1.99pident)
# 17805    14
table(r1.pc.blast.mts1$fscore)

# How many of those have a f-score bigger than 99% ?
sum(r1.pc.blast.mts1$fscore >= 99)
# 17536

# Total number of original sequences
nrow(r1.pc.blast.mts1)
# 22975

# Identify sequence with a lower fscore because uncomplete alignement
others2 <- r1.pc.blast.mts1.99pident[r1.pc.blast.mts1.99pident$fscore < 99,]

sum(others2$send == 1)
# 71

sum (others2$sstart == others2$slen)
# 35

# Roughly the same, only number of seq with f-score >= 99 are different.



test <- "AACAACGCACTATCTGTTCATAGTAGTAGTTGCTGTGAATAGCACTCTGCTGACAATTAATGCAGGAGACTACATCTTCTACACTGACTGGGCTTGGACGTCGTTTGTAGTATTTTCGATATCACAATCGACAATGCTTGCAGTAGGAGCAATTTACTACATGCTCTTTACAGGAGTTCCAGGTACAGCGACATATTACGCTGCTGGTAAGAGATCCCCTGGAAGTAGCATTCAAGTATCCTCGACCGACATTGCCACCTTACATGACGCCTATAGAACCTCAGGTCGGTAAGTTCTACAACAGTCCCGTAGCCTTGGGAGCTGGCGCTGGGGCAGTGCTAACTGTGCCTATTGCAGCGTTAGGTGCGAAACTAAACACG"

nchar(test)
nchar(substring(test, 1, 200))
nchar(substring(test, 201, 380))
