### Taxonomic annotation of amoA probe capture sequences and comparison with amplicon sequencing ###

# Based on sibtran 15_usearch_v5 pipeline #
# Written by alexandre.bagnoud@gmail.com
# version3 09-07-2019

cd ~/OneDrive/Work/Collaboration/henri_probe_capture/3-metagenome

# Raw sequences files given by Henir (in '0-raw_files/')



### Databases

# Alves et al., 2018, Nature Communications, Unifying the global phylogeny and environmental distribution of ammonia-oxidising archaea based on amoA genes

# Set path to database files
db_seq="/Users/Alex/OneDrive/Work/UniVienna/ricardo_amoA_classification/Supplementary_Data_1_Databases_Trees/d_AamoA.db_nr_aln.fasta"
chimera_db="/Users/Alex/OneDrive/Work/UniVienna/ricardo_amoA_classification/Supplementary_Data_1_Databases_Trees/j_AamoA_chimera.ref.db_aln.trim.fasta"
qiime_tax="/Users/Alex/OneDrive/Work/UniVienna/ricardo_amoA_classification/Supplementary_Data_1_Databases_Trees/e_AamoA.db_nr_aln_taxonomy_qiime.txt"
mothur_tax="/Users/Alex/OneDrive/Work/UniVienna/ricardo_amoA_classification/Supplementary_Data_1_Databases_Trees/f_AamoA.db_nr_aln_taxonomy_mothur.txt"



### Programs used

# usearch v8
# qiime1
# mothur

### 1) Dereplicate probe capture reads

usearch_v8.1.1861 --derep_fulllength 0-raw_files/TamoA_H2_17N_hmm_0-1.fasta \
--fastaout 1-derep_mg.fasta --sizeout
# 3 seqs, 3 uniques, 3 singletons (100.0%)

### 2)Cluster OTUs at 97%

usearch_v8.1.1861 -cluster_otus 1-derep_mg.fasta -minsize 1 \
-otus 2-otus_mg.fasta -relabel OTU -sizeout
# 2 OTUs, 0 chimeras (0.0%)


# Select OTU by size (set the min_size variable)
min_size="1"
usearch_v8.1.1861 -sortbysize 2-otus_mg.fasta -fastaout 3-otus_mg_${min_size}.fasta -minsize $min_size
# Sorting 2 sequences

### 4) Discard all sequences that share less than 55% identity with any reference sequences

usearch_v8.1.1861 -usearch_global 3-otus_mg_${min_size}.fasta -db $db_seq -id 0.55 -strand both \
-uc 4a1-uclust_report_47_hmm.txt -matched 4a-otus_mg_min${min_size}_match.fasta -notmatched 4b-otus_mg_min${min_size}_nomatch.fasta
# 100% matched

### 5) UCHIME chimera filtration (using parameters defined by Alves et al., 2018)

# Reverse-complement of chimeric database
usearch_v8.1.1861 -fastx_revcomp $chimera_db -label_suffix _RC -fastaout 5-chimera_db_rc.fasta

usearch_v8.1.1861 -uchime_ref 4a-otus_mg_min${min_size}_match.fasta -db $chimera_db \
-nonchimeras 5a-otus_mg_min${min_size}_nochim_temp.fasta -strand plus -mindiv 1.7 -minh 0.1
# Found 0/2 chimeras (0.0%), 2 not classified (100.0%)

usearch_v8.1.1861 -uchime_ref 5a-otus_mg_min${min_size}_nochim_temp.fasta -db 5-chimera_db_rc.fasta \
-nonchimeras 5b-otus_mg_min${min_size}_nochim.fasta -strand plus -mindiv 1.7 -minh 0.1
# Found 0/2 chimeras (0.0%), 1 not classified (50.0%)

### 6) Annotation of the OTUs

# 6.1) With UCLUST implemented in QIIME1

source activate qiime1
assign_taxonomy.py -i 5b-otus_mg_min${min_size}_nochim.fasta -t $qiime_tax -r $db_seq --similarity 0.7 \
-o 6a-otus_mg_min${min_size}_match_uchimed_uclust_annotation/
source deactivate qiime1

# Number of uniques annotations
less 6a-otus_mg_min${min_size}_match_uchimed_uclust_annotation/5b-otus_mg_min${min_size}_nochim_tax_assignments.txt | cut -f2 | sort -u | wc -l
# 2

# Number of Unassigned OTUs
grep -c Unassigned 6a-otus_mg_min${min_size}_match_uchimed_uclust_annotation/5b-otus_mg_min${min_size}_nochim_tax_assignments.txt
# 0

# Note 1: if you get a lot of unassigned, try to lower the similarity threshold, i.e. 0.8
# Note 2: if the number of Unassigned OTUs is still too big, increase the '-minsize' option in point 6)
# when defining OTUs and repeat the following steps


# 6.2) With Mothur

mothur_dir="6b-otus_mg_min${min_size}_match_uchimed_mother_annotation"

mkdir $mothur_dir
cd $mothur_dir
cp ../5b-otus_mg_min${min_size}_nochim.fasta otus_mg_min${min_size}_match_uchimed.fasta # mothur doesn't like file names starting with "[0-9]-"
mothur "#classify.seqs(fasta=otus_mg_min${min_size}_match_uchimed.fasta, reference=$db_seq, taxonomy=$mothur_tax, cutoff=90)"
cd ..

# Number of uniques annotations
less 6b-otus_mg_min${min_size}_match_uchimed_mother_annotation/otus_mg_min${min_size}_match_uchimed.db_nr_aln_taxonomy_mothur.wang.taxonomy | cut -d ";" -f12 | perl -pe 's/_unclassified\([0-9]*\)//g'| sort -u | wc -l
# 2

# Number of Unassigned OTUs
grep -c "\tunknown;" 6b-otus_mg_min${min_size}_match_uchimed_mother_annotation/otus_mg_min${min_size}_match_uchimed.db_nr_aln_taxonomy_mothur.wang.taxonomy 
# 0

# Note 1: if you get a lot of Unassigned in the annotation, try to lower the similarity cutoff, i.e. 80
# Note 2: if the number of Unassigned OTUs is still too big, increase the '-minsize' option in point 6)
# when defining OTUs and repeat the following steps

### 7) Blast check of chimeras

# Make blast database 
makeblastdb -in $db_seq -out temp5-amoa_db -dbtype nucl

blastn -query 5b-otus_mg_min${min_size}_nochim.fasta -db temp5-amoa_db -outfmt 0 -out 9-blast_report_otus_mg_min${min_size}.txt


### 8) Build an OTU table

# Re-label sequences headers and concatenate probe capture sequences

# Build the OTU table
usearch_v8.1.1861 -usearch_global 0-raw_files/TamoA_H2_17N_hmm_0-1.fasta \
-db 5b-otus_mg_min${min_size}_nochim.fasta -strand plus -id 0.97 -otutabout 10-otu_mg_tab_min${min_size}.txt -matched temp10-raw_reads_mg_match_min${min_size}.fasta
# 3 / 3 mapped to OTUs (100.0%)   













