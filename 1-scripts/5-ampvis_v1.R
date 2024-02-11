###### AMPVIS2

library("ampvis2")
setwd("~/OneDrive/Work/Collaboration/henri_probe_capture/6-ampvis/")

#### 1) Merge OTU tables

### 1.1) Import amplicon data

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

### 1.2) Merge file

amp <- merge(amp.tax, amp.otu)

### 1.3) Aggregate file

amp.agg <- aggregate(amp[,3:5], by = list(amp$tax), FUN = sum)

colnames(amp.agg) <- c("tax", "ampH2_16", "ampH2_17", "ampH2_18")

### 1.4) Import probe capture data

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

### 1.6) Merge file

pc <- merge(pc.tax, pc.otu)

### 1.7) Aggregate file

pc.agg <- aggregate(pc[,3:5], by = list(pc$tax), FUN = sum)

colnames(pc.agg) <- c("tax", "pcH2_16", "pcH2_17", "pcH2_18")


### 1.8) Merge the 3 dataframes

m <- merge(amp.agg, pc.agg, by = "tax", all = TRUE)

m[is.na(m)] <- 0

### 1.9) Save to merge OTU table

write.table(m, "1-merge_OTU_tab.txt", sep = "\t", quote = FALSE, row.names = FALSE)


#### 2) Import files to AMPVIS2


### 2.1) load data

## merged OTU table

otu.tab <- m[,-1]
rownames(otu.tab) <- m$tax
otu.tab$Kingdom <- NA
otu.tab$Phylum <- NA
otu.tab$Class <- NA
otu.tab$Order <- NA
otu.tab$Family <- NA
otu.tab$Genus <- NA
otu.tab$Species <- NA

## Remove unkown OTUs
otu.tab <- otu.tab[rownames(otu.tab) != "unknown",]

## metadata

md <- data.frame(Sample = colnames(otu.tab)[1:6], Source = c("amplicon", "amplicon", "amplicon",
                                                              "probe.capture", "probe.capture", "probe.capture"))

## Import data to AMPVIS2
ampvis <- amp_load(otutable = otu.tab, metadata = md)

#### 3) Plot rarefaction curve

plot1 <- amp_rarecurve(ampvis, color_by = "Source", stepsize = 100) +
    scale_color_manual(name = "Sequencing\napproach", 
                       labels = c("Amplicon seq", "Probe capture"), 
                       values = c("#377eb8", "#e41a1c")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))

print(plot1)

pdf("2-rarefaction_plot.pdf",width=6,height=4) 
print(plot1) 
dev.off()

#### 4) Venn Diagram

library("VennDiagram")

# Generate 2 sets of OTUs

otu.tab$amp <- rowSums(otu.tab[,1:3])
otu.tab$pc <- rowSums(otu.tab[,4:6])

amp.set <- rownames(otu.tab)[otu.tab$amp > 0]
pc.set <- rownames(otu.tab)[otu.tab$pc > 0]


# Chart
plot2 <- venn.diagram(
    x = list(amp.set, pc.set),
    category.names = c("Amplicon\nsequecing" , "Probe capture"),
    filename = NULL,
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = c("#377eb8", "#e41a1c"),
    
    # Numbers
    cex = 1.2,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 1.2,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "sans"
)

ggsave(plot2, file="3-venn_diagram.pdf", device = "pdf", height = 4, width = 4)
