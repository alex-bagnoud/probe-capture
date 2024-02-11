### Taxonomic annotation of amoA probe capture sequences and comparison with amplicon sequencing ###

# Based on sibtran 15_usearch_v5 pipeline #
# Written by alexandre.bagnoud@gmail.com
# version3 09-07-2019

# Raw sequences files (in '0-raw_files/')

# TamoA_47_R1_hmm.fasta = H2-16 (probe capture)
# TamoA_48_R1_hmm.fasta = H2-17 (probe capture)
# TamoA_49_R1_hmm.fasta = H2-18 (probe capture)
# TamoA_H2_17N_hmm.fasta = H2-17 (metagnomic assembly)



### Databases

# Alves et al., 2018, Nature Communications, Unifying the global phylogeny and environmental distribution of ammonia-oxidising archaea based on amoA genes

# Set path to database files
db_seq="0-databases/d_AamoA.db_nr_aln.fasta"
chimera_db="0-databases/j_AamoA_chimera.ref.db_aln.trim.fasta"
qiime_tax="0-databases/e_AamoA.db_nr_aln_taxonomy_qiime.txt"
mothur_tax="0-databases/f_AamoA.db_nr_aln_taxonomy_mothur.txt"



### Programs used

# usearch v8
# qiime1
# mothur

### 1) Concatenate probe capture sequences files
cat 0-raw_files/TamoA_47_R1_hmm.fasta 0-raw_files/TamoA_48_R1_hmm.fasta \
0-raw_files/TamoA_49_R1_hmm.fasta > 0-raw_files/1-tamoA_pc_seq.fasta

### 2) Dereplicate probe capture reads

usearch_v8.1.1861 --derep_fulllength 0-raw_files/1-tamoA_pc_seq.fasta \
--fastaout 1-derep_pc.fasta --sizeout
# 20849 seqs, 14269 uniques, 10324 singletons (72.4%) 

### 3)Cluster OTUs at 97%

usearch_v8.1.1861 -cluster_otus 1-derep_pc.fasta -minsize 1 \
-otus 2-otus_pc.fasta -relabel OTU -sizeout
# 100.0% 238 OTUs, 9 chimeras (0.1%)


# Select OTU by size (set the min_size variable)
min_size="2"
usearch_v8.1.1861 -sortbysize 2-otus_pc.fasta -fastaout 3-otus_pc_${min_size}.fasta -minsize $min_size
# Sorting 198 sequences

### 4) Discard all sequences that share less than 55% identity with any reference sequences

usearch_v8.1.1861 -usearch_global 3-otus_pc_${min_size}.fasta -db $db_seq -id 0.55 -strand both \
-uc 4a1-uclust_report_47_hmm.txt -matched 4a-otus_pc_min${min_size}_match.fasta -notmatched 4b-otus_pc_min${min_size}_nomatch.fasta
# 94.4% matched

### 5) UCHIME chimera filtration (using parameters defined by Alves et al., 2018)

# Reverse-complement of chimeric database
usearch_v8.1.1861 -fastx_revcomp $chimera_db -label_suffix _RC -fastaout 5-chimera_db_rc.fasta

usearch_v8.1.1861 -uchime_ref 4a-otus_pc_min${min_size}_match.fasta -db $chimera_db \
-nonchimeras 5a-otus_pc_min${min_size}_nochim_temp.fasta -strand plus -mindiv 1.7 -minh 0.1
# Found 0/187 chimeras (0.0%), 187 not classified (100.0%)

usearch_v8.1.1861 -uchime_ref 5a-otus_pc_min${min_size}_nochim_temp.fasta -db 5-chimera_db_rc.fasta \
-nonchimeras 5b-otus_pc_min${min_size}_nochim.fasta -strand plus -mindiv 1.7 -minh 0.1
# Found 1/187 chimeras (0.5%), 82 not classified (43.9%)

### 6) Annotation of the OTUs

# 6.1) With UCLUST implemented in QIIME1

source activate qiime1
assign_taxonomy.py -i 5b-otus_pc_min${min_size}_nochim.fasta -t $qiime_tax -r $db_seq --similarity 0.7 \
-o 6a-otus_pc_min${min_size}_match_uchimed_uclust_annotation/
source deactivate qiime1

# Number of uniques annotations
less 6a-otus_pc_min${min_size}_match_uchimed_uclust_annotation/5b-otus_pc_min2_nochim_tax_assignments.txt | cut -f2 | sort -u | wc -l
# 50

# Number of Unassigned OTUs
grep -c Unassigned 6a-otus_pc_min${min_size}_match_uchimed_uclust_annotation/5b-otus_pc_min2_nochim_tax_assignments.txt
# 86

# Note 1: if you get a lot of unassigned, try to lower the similarity threshold, i.e. 0.8
# Note 2: if the number of Unassigned OTUs is still too big, increase the '-minsize' option in point 6)
# when defining OTUs and repeat the following steps


# 6.2) With Mothur

mothur_dir="6b-otus_pc_min${min_size}_match_uchimed_mother_annotation"

mkdir $mothur_dir
cd $mothur_dir
cp ../5b-otus_pc_min${min_size}_nochim.fasta otus_pc_min${min_size}_match_uchimed.fasta # mothur doesn't like file names starting with "[0-9]-"
mothur "#classify.seqs(fasta=otus_pc_min${min_size}_match_uchimed.fasta, reference=$db_seq, taxonomy=$mothur_tax, cutoff=90)"
cd ..

# Number of uniques annotations
less 6b-otus_pc_min${min_size}_match_uchimed_mother_annotation/otus_pc_min${min_size}_match_uchimed.db_nr_aln_taxonomy_mothur.wang.taxonomy | cut -d ";" -f12 | perl -pe 's/_unclassified\([0-9]*\)//g'| sort -u | wc -l
# 84

# Number of Unassigned OTUs
grep -c "\tunknown;" 6b-otus_pc_min${min_size}_match_uchimed_mother_annotation/otus_pc_min${min_size}_match_uchimed.db_nr_aln_taxonomy_mothur.wang.taxonomy 
# 16

# Note 1: if you get a lot of Unassigned in the annotation, try to lower the similarity cutoff, i.e. 80
# Note 2: if the number of Unassigned OTUs is still too big, increase the '-minsize' option in point 6)
# when defining OTUs and repeat the following steps

### 7) Check the chimeric nature of the unassigned OTUs.

# Make blast database 
makeblastdb -in $db_seq -out temp5-amoa_db -dbtype nucl

# Make the fasta flat (for easier selection of OTUs)
less 5b-otus_pc_min${min_size}_nochim.fasta | \
awk -v RS='>'     -v FS="\n"     -v OFS=""     -v ORS="" '{ if (NR > 1) { printf ">%s\n",$1; $1=""; printf "%s\n",$0 } }' \
> 5b-otus_pc_min${min_size}_nochim_flat.fasta

# List of unassigned OTUs from UCLUST annotation
grep Unassigned 6a-otus_pc_min${min_size}_match_uchimed_uclust_annotation/5b-otus_pc_min2_nochim_tax_assignments.txt | cut -f1 > temp6a-unassigned_uclust_otu_list.txt

# List of unassigned OTUs from MOTHUR annotation
grep "\tunknown;" 6b-otus_pc_min${min_size}_match_uchimed_mother_annotation/otus_pc_min${min_size}_match_uchimed.db_nr_aln_taxonomy_mothur.wang.taxonomy | cut -f1 > temp6b-unassigned_mothur_otu_list.txt

# Extract the unassigned sequences from the flat fast file for UCLUST annotations
grep -A1 -f temp6a-unassigned_uclust_otu_list.txt 5b-otus_pc_min${min_size}_nochim_flat.fasta \
| grep -v "^--$" > 7a_unassigned_pc_uclust_min${min_size}_nochim_amoa_otus.fasta

# Extract the unassigned sequences from the flat fast file for MOTHUR annotations
grep -A1 -f temp6b-unassigned_mothur_otu_list.txt 5b-otus_pc_min${min_size}_nochim_flat.fasta \
| grep -v "^--$" > 7b_unassigned_pc_mothur_min${min_size}_nochim_amoa_otus.fasta

# Run blast (UCLUST)
blastn -query 7a_unassigned_pc_uclust_min${min_size}_nochim_amoa_otus.fasta -db temp5-amoa_db -outfmt 0 -out 8a-blast_report_unassigned_uclust_min${min_size}_otus.txt

# Run blast (MOTHUR)
blastn -query 7b_unassigned_pc_mothur_min${min_size}_nochim_amoa_otus.fasta -db temp5-amoa_db -outfmt 0 -out 8b-blast_report_unassigned_mothur_min${min_size}_otus.txt

# Carefully review the BLAST alignement and determine if the OTUs sequences are artefacts or not.
# If some are, manually remove them from the file '6b-otus_min${min_size}_match_uchimed_flat.fasta'
# and save it as '10_clean_otus_min${min_size}.fasta'

# Here a lot of the unassigned sequences from uclust have multiple perfect/very good matches with Ricardo's db.
# Not sure what causes this...
# This is also true for some unassigned OTUs from Mothur
# Maybe these sequences are very conserved and thus hard to classify...


# Quick check the blast alignment of all OTUs
blastn -query 5b-otus_pc_min${min_size}_nochim.fasta -db temp5-amoa_db -outfmt 0 -out 9a-blast_report_otus_pc_min${min_size}.txt
blastn -query 5b-otus_pc_min${min_size}_nochim.fasta -db temp5-amoa_db -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send slen evalue bitscore" \
-out 9b-blast_report_otus_pc_min${min_size}_tab.txt



### 8) Build an OTU table

# Re-label sequences headers and concatenate probe capture sequences
less 0-raw_files/TamoA_47_R1_hmm.fasta | perl -pe 's/^>.*$/>H2_16/g' > 0-raw_files/TamoA_47_R1_hmm_sample_labels.fasta
less 0-raw_files/TamoA_48_R1_hmm.fasta | perl -pe 's/^>.*$/>H2_17/g' > 0-raw_files/TamoA_48_R1_hmm_sample_labels.fasta
less 0-raw_files/TamoA_49_R1_hmm.fasta | perl -pe 's/^>.*$/>H2_18/g' > 0-raw_files/TamoA_49_R1_hmm_sample_labels.fasta

cat 0-raw_files/TamoA_47_R1_hmm_sample_labels.fasta 0-raw_files/TamoA_48_R1_hmm_sample_labels.fasta \
0-raw_files/TamoA_49_R1_hmm_sample_labels.fasta > 0-raw_files/2-tamoA_pc_seq_sample_labels.fasta

# Build the OTU table
usearch_v8.1.1861 -usearch_global 0-raw_files/2-tamoA_pc_seq_sample_labels.fasta \
-db 5b-otus_pc_min${min_size}_nochim.fasta -strand plus -id 0.97 -otutabout 10-otu_pc_tab_min${min_size}.txt -matched temp10-raw_reads_pc_match_min${min_size}.fasta
# 20688 / 20849 mapped to OTUs (99.2%) 














