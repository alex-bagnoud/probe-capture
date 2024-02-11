#!/usr/bin/env Rscript

# Author: alexandre.bagnoud@gmail.com (April 2018)

# This script fuses fasta sequences from two paired fasta files, for all paired fasta file present in one directory
# Forward files should be named 'sample_R1.fasta' and reverse files 'sample_R2rc.fasta' (rc stands for reverse-complement).
# The script fuses reads based on a similar sequences header, but all character after the first space are ignored.    
# The script also replaces the header by the sample name, the one contained in the files names.

library("optparse")

option_list = list(
    make_option(c("-i", "--folderIN"), type="character", default="2-qc_filtered_fasta_files/", 
                help="Input folder that contains all fasta files to fuse. Each sample should be splitted into one '_R1.fasta' and another '_R2rc.fasta'.", metavar="character"),
    make_option(c("-o", "--folderOUT"), type="character", default="3-fused_fasta", 
                help="Output folder", metavar="character")

); 

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Libraries
library("Biostrings")
library("stringr")

# Functions
fastaToDf <- function(fastaFile){
    dnaSeq <- readBStringSet(fastaFile)
    fasta_df <- data.frame(header = names(dnaSeq), sequence = paste(dnaSeq))
}

dfToFasta <- function(header, seq, file){
    sequence = BStringSet(seq)
    names(sequence) <- header
    writeXStringSet(sequence, file)
}


# Variables
path <- opt$folderIN
output.dir <- opt$folderOUT

dir.create(output.dir)
file_list <- sort(list.files(path, pattern="_R1.fasta", full.names = TRUE))
file_list <- gsub("//", "/", file_list)

for (file in file_list){
    r1_file_name <- file
    r2_file_name <- sub("R1", "R2rc", r1_file_name)
    sample_name <- str_split_fixed(basename(r1_file_name), "_", 2)[1,1]
    output.file <- paste0(output.dir, "/", sample_name, "_fused.fasta")
    
    if (!file.exists(r2_file_name)){
        file.create(output.file)
        
    } else {
        
        r1 <- fastaToDf(r1_file_name)
        r2 <- fastaToDf(r2_file_name)
        
        r1$header <- str_split_fixed(r1$header, " ", 2)[,1]
        r2$header <- str_split_fixed(r2$header, " ", 2)[,1]
        
        r1r2 <- merge(r1, r2, by.x = "header", by.y = "header")
        
        r1r2$fused <- paste0(r1r2$sequence.x, r1r2$sequence.y)
        r1r2$sample <- rep(sample_name, nrow(r1r2))
        
        dfToFasta(r1r2$sample, r1r2$fused, output.file)
    }

}
