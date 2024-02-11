#### Usearch pipeline for analysing non-overlapping paired-ends amplicons of archaeal amoA ####

# Written by alexandre.bagnoud@gmail.com
# version3 11-07-2019
# based on sibtran pipeline (15_usearch_v5)


# NOTE: The pipeline proposes two softwares to annotate the OTUs, UCLUST (implemented in QIIME1)
# and MOTHUR. The first one is better as detected artefact sequences, by annotating them
# as 'Unassigned'. MOTHUR doesn't point them out, as it give them a well-defined annotation.
# For this reason, I would recommend to use UCLUST to annotate OTUs.



### Programs used

# cutadapt version 1.13
# usearch v8
# R version (and packages)
# mothur
# qiime1




### Databases

# Alves et al., 2018, Nature Communications, Unifying the global phylogeny and environmental distribution of ammonia-oxidising archaea based on amoA genes

# Set path to database files
db_seq="0-databases/d_AamoA.db_nr_aln.fasta"
qiime_tax="0-databases/e_AamoA.db_nr_aln_taxonomy_qiime.txt"
mothur_tax="0-databases/f_AamoA.db_nr_aln_taxonomy_mothur.txt"
chimera_db="0-databases/j_AamoA_chimera.ref.db_aln.trim.fasta"



### Other scripts

# fuse_fasta_paired_reads.R (executable from bash)
# otu_plots_and_tables.R (to be run in R)



### 1) Removing non-biological sequences from fastq files
# Here only the first 'G' still belongs to the primer in the forward primer

## set variables:
fastq_path="0-raw_files/"
r1="*R1.fastq"
r2="*R2.fastq"
trim_dir="temp1-trimmed_fastq_files/"


mkdir $trim_dir

## Forward read: remove first 'G'
for file in $fastq_path$r1
do
sample="${file##*/}"
cutadapt -g G -o $trim_dir$sample $file
done

cat ${trim_dir}*${r1} > ${trim_dir}all_R1.fastq

## Reverse reads: simply copy files
cp $fastq_path$r2 $trim_dir

cat ${trim_dir}*${r2} > ${trim_dir}all_R2.fastq



### 2) Quality filtration of reads

# Which parameters quality to choose for R1 reads?
# The ee-cutoff of 500 will display the proportion of reads regardless their quality in the different length bins
usearch_v8.1.1861 -fastq_eestats2 ${trim_dir}all_R1.fastq -length_cutoffs 150,300,10 -ee_cutoffs 0.5,1,2,3,4,5,500 \
-output 1a-eestats_r1.txt
cat 1a-eestats_r1.txt

# Which parameters quality to choose for R2 reads?
# The ee-cutoff of 500 will display the proportion of reads regardless their quality in the different length bins
usearch_v8.1.1861 -fastq_eestats2 ${trim_dir}all_R2.fastq -length_cutoffs 150,300,10 -ee_cutoffs 0.5,1,2,3,4,5,500 \
-output 1b-eestats_r2.txt
cat 1b-eestats_r2.txt

# Note: you can be stringent here, as there reads will only be used for defining OTUs.
# Once OTUs are defined, all the raw reads will be truncated at the same length and fused,
# and mapped to the OTUs in order to compute an OTU table.


# Quality filtration
qc_filt_dir="temp2-qc_filtered_fasta_files/"
mkdir $qc_filt_dir

# Quality parameters
r1_length="200"
r2_length="180"
r1_maxee="0.5"
r2_maxee="0.5"

# R1 reads
for file in $fastq_path$r1
do
sample="$(basename "$file" .fastq)"
usearch_v8.1.1861 -fastq_filter $trim_dir${sample}.fastq --fastq_trunclen $r1_length --fastq_maxee $r1_maxee --fastaout $qc_filt_dir${sample}.fasta
done

# R2 reads
for file in $fastq_path$r2
do
sample="$(basename "$file" .fastq)"
usearch_v8.1.1861 -fastq_filter $trim_dir${sample}.fastq --fastq_trunclen $r2_length --fastq_maxee $r2_maxee --fastaout $qc_filt_dir${sample}.fasta
done



### 3) Reverse complement of R2 sequences

for file in ${qc_filt_dir}*R2.fasta
do
sample="$(basename "$file" .fasta)"
echo $sample
usearch_v8.1.1861 --fastx_revcomp $file --fastaout $qc_filt_dir${sample}rc.fasta
done



### 4) Fuse reads and rename reads with their sample names

./fuse_fasta_paired_reads.R -i $qc_filt_dir -o "temp3-fused_fasta"

# Note: some fused fasta files might be empty due to the fact that some samples with a very low sequencing
# depth might loose all their sequences. This sample might be discarded for the OTU selection,
# or less stringent quality parameters might be chosen. In this case the following command
# might be useful
# "usearch -fastq_eestats2 temp1-trimmed_fastq_files/sample.fastq -length_cutoffs 150,300,10 -ee_cutoffs 0.5,1,2,3,4,5,500"

# Concatenate all these files into one
cat "temp3-fused_fasta"/*.fasta > temp4-filtered_merged_reads.fasta

# Tracking reads along this pipeline
echo -e "sample\traw\tqcR1\tqcR2\tfused" > 2-track.txt
for file in temp3-fused_fasta/*
do 
sample=$(basename "$file" _fused.fasta)
count1=$(wc -l $trim_dir${sample}_R1.fastq | tr -s " " | cut -f2 -d " ")
count1b=`echo "$count1/4" | bc `
count2=$(grep -c "^>" $qc_filt_dir${sample}_R1.fasta)
count3=$(grep -c "^>" $qc_filt_dir${sample}_R2.fasta)
count4=$(grep -c "^>" temp3-fused_fasta/${sample}_fused.fasta)
echo -e "${sample}\t${count1b}\t${count2}\t${count3}\t${count4}" >> 2-track.txt
done
cat 2-track.txt

# sample	raw	qcR1	qcR2	fused
# H2-16	58006	43719	44308	34995
# H2-17	27565	20274	20462	15677
# H2-18	35955	27098	27135	21379


### 5) Dereplicate reads

usearch_v8.1.1861 --derep_fulllength temp4-filtered_merged_reads.fasta \
--fastaout 3-filterd_merged_derep.fasta --sizeout
# 72051 seqs, 41452 uniques, 37157 singletons (89.6%)


### 6) Cluster OTUs at 97%

usearch_v8.1.1861 -cluster_otus 3-filterd_merged_derep.fasta -minsize 1 \
-otus 4a-otus.fasta -relabel OTU -sizeout
# 201 OTUs, 6304 chimeras (15.2%)

# Select OTU by size (set the min_size variable)
min_size="2"
usearch_v8.1.1861 -sortbysize 4a-otus.fasta -fastaout 4b-otus_min${min_size}.fasta -minsize $min_size
# Sorting 119 sequences


### 7) Discard all sequences that share less than 55% identity with any reference sequences

usearch_v8.1.1861 -usearch_global 4b-otus_min${min_size}.fasta -db $db_seq -id 0.55 -strand plus \
-uc 5a-uclust_report.txt -matched 5b-otus_min${min_size}_match.fasta -notmatched 5c-otus_min${min_size}_nomatch.fasta
# 100.0% matched



### 8) UCHIME chimera filtration (using parameters defined by Alves et al., 2018)

usearch_v8.1.1861 -uchime_ref 5b-otus_min${min_size}_match.fasta -db $chimera_db \
-nonchimeras 6a-otus_min${min_size}_match_uchimed.fasta -strand plus -mindiv 1.7 -minh 0.1
# Found 15/119 chimeras (12.6%), 52 not classified (43.7%)



### 9) Annotation of the OTUs

# 9.1) With UCLUST implemented in QIIME1

source activate qiime1
assign_taxonomy.py -i 6a-otus_min${min_size}_match_uchimed.fasta -t $qiime_tax -r $db_seq --similarity 0.9 \
-o 7a-otus_min${min_size}_match_uchimed_uclust_annotation/
source deactivate qiime1

# Number of uniques annotations
less 7a-otus_min${min_size}_match_uchimed_uclust_annotation/6a-otus_min${min_size}_match_uchimed_tax_assignments.txt | cut -f2 | sort -u | wc -l
# 61

# Number of Unassigned OTUs
grep -c Unassigned 7a-otus_min${min_size}_match_uchimed_uclust_annotation/6a-otus_min${min_size}_match_uchimed_tax_assignments.txt
# 1

# Note 1: if you get a lot of unassigned, try to lower the similarity threshold, i.e. 0.8
# Note 2: if the number of Unassigned OTUs is still too big, increase the '-minsize' option in point 6)
# when defining OTUs and repeat the following steps

# 9.2) With Mothur

mothur_dir="7b-otus_min${min_size}_match_uchimed_mothur_annotation"

mkdir $mothur_dir
cd $mothur_dir
cp ../6a-otus_min${min_size}_match_uchimed.fasta otus_min${min_size}_match_uchimed.fasta # mothur doesn't like file names starting with "[0-9]-"
mothur "#classify.seqs(fasta=otus_min${min_size}_match_uchimed.fasta, reference=$db_seq, taxonomy=$mothur_tax, cutoff=90)"
cd ..

# Note 1: if you get a lot of Unassigned in the annotation, try to lower the similarity cutoff, i.e. 80
# Note 2: if the number of Unassigned OTUs is still too big, increase the '-minsize' option in point 6)
# when defining OTUs and repeat the following steps

# Number of uniques annotations
less 7b-otus_min${min_size}_match_uchimed_mothur_annotation/otus_min${min_size}_match_uchimed.db_nr_aln_taxonomy_mothur.wang.taxonomy | cut -d ";" -f12 | perl -pe 's/_unclassified\([0-9]*\)//g'| sort -u | wc -l
# 68

# Number of Unassigned OTUs
grep -c "\tunknown;" 7b-otus_min${min_size}_match_uchimed_mothur_annotation/otus_min${min_size}_match_uchimed.db_nr_aln_taxonomy_mothur.wang.taxonomy
# 0


### 10) Check the chimeric nature of the unassigned OTUs.

# Make blast database 
makeblastdb -in $db_seq -out temp5-amoa_db -dbtype nucl

# Make the fasta flat (for easier selection of OTUs)
less 6a-otus_min${min_size}_match_uchimed.fasta | \
awk -v RS='>'     -v FS="\n"     -v OFS=""     -v ORS="" '{ if (NR > 1) { printf ">%s\n",$1; $1=""; printf "%s\n",$0 } }' \
> 6b-otus_min${min_size}_match_uchimed_flat.fasta

# List of unassigned OTUs from UCLUST annotation
grep Unassigned 7a-otus_min${min_size}_match_uchimed_uclust_annotation/6a-otus_min${min_size}_match_uchimed_tax_assignments.txt | cut -f1 > temp6a-unassigned_uclust_otu_list.txt

# List of unassigned OTUs from MOTHUR annotation
grep "\tunknown;" 7b-otus_min${min_size}_match_uchimed_mothur_annotation/otus_min${min_size}_match_uchimed.db_nr_aln_taxonomy_mothur.wang.taxonomy | cut -f1 > temp6b-unassigned_mothur_otu_list.txt

# Extract the unassigned sequences from the flat fast file for UCLUST annotations
grep -A1 -f temp6a-unassigned_uclust_otu_list.txt 6b-otus_min${min_size}_match_uchimed_flat.fasta \
| grep -v "^--$" > 8a_unassigned_uclust_min${min_size}_nochim_amoa_otus.fasta

# Extract the unassigned sequences from the flat fast file for MOTHUR annotations
grep -A1 -f temp6b-unassigned_mothur_otu_list.txt 6b-otus_min${min_size}_match_uchimed_flat.fasta \
| grep -v "^--$" > 8b_unassigned_mothur_min${min_size}_nochim_amoa_otus.fasta

# Run blast (UCLUST)
blastn -query 8a_unassigned_uclust_min${min_size}_nochim_amoa_otus.fasta -db temp5-amoa_db -outfmt 0 -out 9a-blast_report_unassigned_uclust_min${min_size}_otus.txt

# Run blast (MOTHUR)
blastn -query 8b_unassigned_mothur_min${min_size}_nochim_amoa_otus.fasta -db temp5-amoa_db -outfmt 0 -out 9b-blast_report_unassigned_mothur_min${min_size}_otus.txt

# Carefully review the BLAST alignement and determine if the OTUs sequences are artefacts or not.
# If some are, manually remove them from the file '6b-otus_min${min_size}_match_uchimed_flat.fasta'
# and save it as '10_clean_otus_min${min_size}.fasta'

# Quick check the blast alignment of all OTUs
blastn -query 6a-otus_min${min_size}_match_uchimed.fasta -db temp5-amoa_db -outfmt 0 -out 9c-blast_report_otus_min${min_size}.txt

### 10) Build an OTU table

# 10.1) Merge all raw paired-end reads together without any quality trimming (except for read length,
# which should be the same as in point 2)
# These sequences will be later used by usearch to create an OTU table, by mapping all these reads
# to the OTUs.

# Raw reads merged at the same length as qc-filtered reads

trunc_dir="temp7-truncated_reads/"

mkdir $trunc_dir

# Trim R1 reads at the right length
for file in $fastq_path$r1
do
sample="$(basename "$file" .fastq)"
usearch_v8.1.1861 -fastq_filter $trim_dir${sample}.fastq --fastq_trunclen $r1_length --fastaout $trunc_dir${sample}.fasta
done

# Trim R2 reads at the right length
for file in $fastq_path$r2
do
sample="$(basename "$file" .fastq)"
usearch_v8.1.1861 -fastq_filter $trim_dir${sample}.fastq --fastq_trunclen $r2_length --fastaout $trunc_dir${sample}.fasta
done

# Reverse complement of R2 sequences

for file in ${trunc_dir}*R2.fasta
do
sample="$(basename "$file" .fasta)"
usearch_v8.1.1861 --fastx_revcomp $file --fastaout $trunc_dir${sample}rc.fasta
done


#Fuse reads and rename reads with their sample names

./fuse_fasta_paired_reads.R -i $trunc_dir -o "temp8-truncated_fused_fasta"

# Concatenate all these files into one and replace all '-' by '_' ('_' not supperted by usearch)
cat temp8-truncated_fused_fasta/*.fasta | perl -pe 's/\-/_/g' > temp9-merged_truncated_raw_reads.fasta


# 10.2) Build the actual OTU table  

usearch_v8.1.1861 -usearch_global temp9-merged_truncated_raw_reads.fasta \
-db 10-clean_otus_min${min_size}.fasta -strand plus -id 0.97 -otutabout 11-otu_tab_min${min_size}.txt -matched temp10-raw_reads_match_min${min_size}.fasta
# 100310 / 121408 mapped to OTUs (82.6%) 

# Track the reads along the pipeline
echo -e "sample\traw\tqcR1\tqcR2\tfused\totutab" > 12-track_min${min_size}.txt
for file in temp3-fused_fasta/*
do 
sample=$(basename "$file" _fused.fasta)
sample2=$(basename "$file" _fused.fasta | sed -e 's/\-/_/g')
count1=$(wc -l $trim_dir${sample}_R1.fastq | tr -s " " | cut -f2 -d " ")
count1b=`echo "$count1/4" | bc `
count2=$(grep -c "^>" $qc_filt_dir${sample}_R1.fasta)
count3=$(grep -c "^>" $qc_filt_dir${sample}_R2.fasta)
count4=$(grep -c "^>" temp3-fused_fasta/${sample}_fused.fasta)
count5=$(grep -c "^>${sample2}" temp10-raw_reads_match_min${min_size}.fasta)
echo -e "${sample}\t${count1b}\t${count2}\t${count3}\t${count4}\t${count5}" >> 12-track_min${min_size}.txt
done
cat 12-track_min${min_size}.txt

# sample	raw	qcR1	qcR2	fused	otutab
# H2-16	58006	43719	44308	34995	47934
# H2-17	27565	20274	20462	15677	23747
# H2-18	35955	27098	27135	21379	28629

