setwd('~/OneDrive/Work/Collaboration/henri_probe_capture/4-pooling/')

##### 1) Import probe capture data

# tax file
utax.pc <- read.table("../2-probe_capture/6a-otus_pc_min2_match_uchimed_uclust_annotation/5b-otus_pc_min2_nochim_tax_assignments.txt", sep = "\t")
colnames(utax.pc)[1:2] <- c("header", "tax") 
utax.pc$tax <- gsub("^.*;", "", utax.pc$tax)
utax.pc$header <- gsub(";.*;", "", utax.pc$header)

mtax.pc <- read.table("../2-probe_capture/6b-otus_pc_min2_match_uchimed_mother_annotation/otus_pc_min2_match_uchimed.db_nr_aln_taxonomy_mothur.wang.taxonomy", sep = "\t")
colnames(mtax.pc) <- c("header", "tax") 
mtax.pc$tax <- gsub("^.*;", "", (gsub(";$", "", mtax.pc$tax)))
mtax.pc$header <- gsub(";.*;", "", mtax.pc$header)

# otu file
otu.tab.pc <- read.table("../2-probe_capture/10-otu_pc_tab_min2.txt", sep = "\t", header = TRUE, comment.char = "@")
colnames(otu.tab.pc)[1] <- "header"

# sequences

library("Biostrings")

fastaToDf <- function(fastaFile){
    dnaSeq <- readBStringSet(fastaFile)
    fasta_df <- data.frame(header = names(dnaSeq), sequence = paste(dnaSeq))
}

fasta.pc <- fastaToDf("../2-probe_capture/5b-otus_pc_min2_nochim.fasta")
fasta.pc$header <- gsub(";.*;", "", fasta.pc$header)
fasta.pc$sequence <- as.character(fasta.pc$sequence)
fasta.pc$seq.len <- nchar(fasta.pc$sequence)

##### 2) Merge the 4 probe capture files

pc1 <- merge(utax.pc, otu.tab.pc)
pc2 <- merge(pc1, fasta.pc)
pc3 <- merge(pc2, mtax.pc, by = "header")

colnames(pc3)[c(2,10)] <- c("utax", "mtax")

write.table(pc3, "1-tax_comparison.txt", sep = "\t", quote = FALSE, row.names = FALSE)

##### 3) Distribution of read length "unassigned" vs "assigned

pc2$assignement <- "Assigned"

for (r in 1:nrow(pc2)) {
    if (pc2$tax[r] == "Unassigned") {
        pc2$assignement[r] <- "Unassigned"
    }
}

library("reshape2")
library("ggplot2")

pc2.m <- melt(data = pc2, id.vars = c("assignement"), measure.vars = c("seq.len"))

plot1 <- ggplot(pc2.m, aes(x = assignement, y = value)) +
    geom_boxplot() +
    geom_point()

print(plot1)

pdf("2-amoa_length_distribution.pdf", width = 5, height = 5)
print(plot1)
dev.off()


##### 4) Analyse BLAST results to assess taxonomy annotations

blast.pc <- read.table("../2-probe_capture/9b-blast_report_otus_pc_min2_tab.txt", sep = "\t")
colnames(blast.pc) <- c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "slen", "evalue", "bitscore")

tax.db <- read.table("../../../UniVienna/ricardo_amoA_classification/Supplementary_Data_1_Databases_Trees/e_AamoA.db_nr_aln_taxonomy_qiime.txt", sep = "\t")
colnames(tax.db) <- c("accession", "tax.hit")
blast.pc2 <- merge(blast.pc, tax.db, by.x = "saccver", by.y = "accession")

blast.pc2$qlen <- NA
blast.pc2$utax <- NA
blast.pc2$mtax <- NA

for (r in 1:nrow(pc3)) {
    log <- grepl(paste0(pc3$header[r], ";"), blast.pc2$qaccver)
    blast.pc2$qlen[log] <- pc3$seq.len[r]
    blast.pc2$utax[log] <- pc3$utax[r]
    blast.pc2$mtax[log] <- pc3$mtax[r]
}

blast.pc2$qcov <- abs(blast.pc2$qstart - blast.pc2$qend) / blast.pc2$qlen
blast.pc2$scov <- abs(blast.pc2$sstart - blast.pc2$send) / blast.pc2$slen

blast.pc2$tax.hit <- gsub("^.*;", "", blast.pc2$tax.hit)


write.table(blast.pc2, "3-blast_pc_with_tax.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# After having look at this file, it seems that mothur annotation are more correct than uclust ones.
# For the OTU that have only been annotated with mother, annotation seem correct based on best blast hits
# For the subsequent analyses, only mother annotation will be kept.

##### 4) Bubble plots

### 4.1) Import amplicon data

# tax file
amp.tax <- read.table("../1-amplicon_seq/7b-otus_min2_match_uchimed_mothur_annotation/otus_min2_match_uchimed.db_nr_aln_taxonomy_mothur.wang.taxonomy",
                      sep = "\t")
colnames(amp.tax) <- c("header", "tax") 
library("stringr")
col_number <- 10
vector_subset <- str_split_fixed(amp.tax$tax, ";", col_number)[,]
amp.tax$tax <- gsub("_unclassified;", "", gsub("\\(.*\\)", "", vector_subset[,10]))
amp.tax$header <- gsub(";.*;", "", amp.tax$header)

# otu file
amp.otu <- read.table("../1-amplicon_seq/11-otu_tab_min2.txt", sep = "\t", comment.char = "@", header = TRUE)
colnames(amp.otu)[1] <- "header"

### 4.2) Merge file

amp <- merge(amp.tax, amp.otu)

### 4.3) Aggregate file

amp.agg <- aggregate(amp[,3:5], by = list(amp$tax), FUN = sum)

### 4.4) Compute relative abundance

amp.agg.counts <- colSums(amp.agg[,-1])
for (c in 2:ncol(amp.agg)) {
    amp.agg[,c] <- amp.agg[,c] / amp.agg.counts[c-1]
}
# Check that the sum of relative abundances for each sample is 1.
colSums(amp.agg[,-1])

colnames(amp.agg) <- c("tax", "ampH2_16", "ampH2_17", "ampH2_18")

### 4.5) Import probe capture data

# tax file
pc.tax <- read.table("../2-probe_capture/6b-otus_pc_min2_match_uchimed_mother_annotation/otus_pc_min2_match_uchimed.db_nr_aln_taxonomy_mothur.wang.taxonomy",
                      sep = "\t")
colnames(pc.tax) <- c("header", "tax") 
library("stringr")
col_number <- 10
vector_subset <- str_split_fixed(pc.tax$tax, ";", col_number)[,]
pc.tax$tax <- gsub("_unclassified;", "", gsub("\\(.*\\)", "", vector_subset[,10]))
pc.tax$header <- gsub(";.*;", "", pc.tax$header)

# otu file
pc.otu <- read.table("../2-probe_capture/10-otu_pc_tab_min2.txt", sep = "\t", comment.char = "@", header = TRUE)
colnames(pc.otu)[1] <- "header"

### 4.6) Merge file

pc <- merge(pc.tax, pc.otu)

### 4.7) Aggregate file

pc.agg <- aggregate(pc[,3:5], by = list(pc$tax), FUN = sum)

### 4.8) Compute relative abundance

pc.agg.counts <- colSums(pc.agg[,-1])
for (c in 2:ncol(pc.agg)) {
    pc.agg[,c] <- pc.agg[,c] / pc.agg.counts[c-1]
}
# Check that the sum of relative abundances for each sample is 1.
colSums(pc.agg[,-1])

colnames(pc.agg) <- c("tax", "pcH2_16", "pcH2_17", "pcH2_18")


### 4.9) Import probe capture data

# tax file
mg.tax <- read.table("../3-metagenome/6b-otus_mg_min1_match_uchimed_mother_annotation/otus_mg_min1_match_uchimed.db_nr_aln_taxonomy_mothur.wang.taxonomy",
                     sep = "\t")
colnames(mg.tax) <- c("header", "tax") 
library("stringr")
col_number <- 10
vector_subset <- str_split_fixed(mg.tax$tax, ";", col_number)[,]
mg.tax$tax <- gsub("_unclassified;", "", gsub("\\(.*\\)", "", vector_subset[,10]))
mg.tax$header <- gsub(";.*;", "", mg.tax$header)

# otu file
mg.otu <- read.table("../3-metagenome/10-otu_mg_tab_min1.txt", sep = "\t", comment.char = "@", header = TRUE)
colnames(mg.otu)[1] <- "header"

### 4.10) Merge file

mg <- merge(mg.tax, mg.otu)

### 4.11) Aggregate file

mg.agg <- aggregate(mg[,3], by = list(mg$tax), FUN = sum)

### 4.12) Compute relative abundance

mg.agg.counts <- sum(mg.agg[,-1])
for (c in 2:ncol(mg.agg)) {
    mg.agg[,c] <- mg.agg[,c] / mg.agg.counts[c-1]
}
# Check that the sum of relative abundances for each sample is 1.
sum(mg.agg[,-1])

colnames(mg.agg) <- c("tax", "mgH3_16")

### 4.13) Merge the 3 dataframes

m1 <- merge(amp.agg, pc.agg, by = "tax", all = TRUE)
m2 <- merge(m1, mg.agg, by = "tax", all = TRUE)

### 4.14) Bubble plot

m2$col <- gsub("-[0-9].*$", "", m2$tax)
m2$col <- gsub("_OTU[0-9]", "", m2$col)

m2[m2 == 0] <- NA

mm <- melt(m2, id.vars = c("tax", "col"))

bubble_plot <- ggplot(mm,aes(variable,tax)) +
    geom_point(aes(size=value, fill=col),shape=21,color="black") +
    theme(panel.grid.major=element_line(linetype=1,color="grey"),
          axis.text.x=element_text(angle=90,hjust=1,vjust=0),
          panel.background = element_blank()) +
    ylab("Taxonomic bins") +
    xlab("Samples") +
    scale_fill_brewer(palette="Set1", name="Taxonomic\nclade") +
    scale_size(name = "Relative\nabundance")

print(bubble_plot)

pdf("4-bubble_plot.pdf", height = 12, width = 6)
print(bubble_plot)
dev.off() 




