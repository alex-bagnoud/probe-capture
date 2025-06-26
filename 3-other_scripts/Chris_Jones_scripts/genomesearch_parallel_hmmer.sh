#!/bin/bash

usage="$(basename "$0") [-h] [-f <alignment_file>.fasta] [-d <path/to/directory>] [-o <outfile.fasta>]\n\n

Program for searching complete genomes for marker gene sequences using HMMER.\n
Ensure that files to be searched are gzipped .fna files from NCBI (extensions are <whatever_genomes>.fna.gz)\n\n

Inputs:\n

    -f  amino acid alignment for generating HMM profile, in fasta format\n
    -d  target directory containing all genomes\n
    -o	name of output file (fasta of retrieved gene sequences)\n"


while getopts hf:d:o: option
do
	case "${option}"
	in
	h) echo -e $usage; exit 1 ;;
	f) INFILE=${OPTARG};;
	d) DIR=${OPTARG};;
	o) OUTFILE=${OPTARG} ;;
   \?) echo -e $usage; exit 1 ;;
	esac
done
	
#for running the hmmer_ch.sh script in parallel. Output files are moved
#to the respective HMMOUT/NUCOUT directories, then the NUCOUT files are
#concatenated. BE sure to change the model files as needed

echo "Making directories NUCOUT and HMMOUT to store result files"

mkdir ./NUCOUT ./HMMOUT

echo $INFILE $DIR $OUTFILE

echo "Generating HMM profile from $INFILE"

hmmbuild tmp.hmm $INFILE

echo "Searching gzipped fna files in $DIR"

find $DIR -name "*.gz" | parallel ./genomesearch_hmmer.sh -f {} -m ./tmp.hmm

for f in ./NUCOUT/*_nucout.fna; do cat "$f" >> ./tmp_nucout.fna; done

sed "s|$DIR||g" ./tmp_nucout.fna | sed "s/>\/GCA/>GCA/g" | sed "s/:>/:/g" | sed ':a;N;$!ba;s/\n\n/\n/g'  > $OUTFILE

rm -r ./tmp.hmm ./tmp_nucout.fna ./NUCOUT

mv ./HMMOUT ./HMMOUT_$OUTFILE
