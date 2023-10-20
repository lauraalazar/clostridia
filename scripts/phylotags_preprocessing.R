#######################
# Initial commands ----
#######################

# load libraries
library(phyloseq)
library(ape)
library(qiime2R)
library(Biostrings)
library(DESeq2)
library(yarrr)
library(stringr)
library(purrr)
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(taxize)
library(RColorBrewer)
library(microbiome)
library(knitr)
library(ggpubr)
library(cowplot)
library(data.table)
library(broom)
library(gridExtra)
library(car)

# define the working directory
setwd("C:/Documentos Disco D/Vidarium/GitHub/clostridia/data")


###############
# METADATA ----
###############

# Metadata files with info for each library
ssu_meta = read.table("metadata/16S.files", header = F)
dnaK_meta = read.table("metadata/dnaK.files", header = F)
gyrB_meta = read.table("metadata/gyrB.files", header = F)

# Edit and remove library names
ssu_meta$lib <- sapply(strsplit(as.character(ssu_meta[,2]), "\\_"), `[`, 1)
dnaK_meta$lib <- sapply(strsplit(as.character(dnaK_meta[,2]), "\\_"), `[`, 1)
gyrB_meta$lib <- sapply(strsplit(as.character(gyrB_meta[,2]), "\\_"), `[`, 1)

# Metadata BMI
sample_info <- read.table("metadata/bsp_classification.txt",header=T)
rownames(sample_info) <- sapply(strsplit(as.character(sample_info$ID),"_H"), `[`, 1)

# Metadata for design variables
microbio_selected <- read.table("metadata/microbio_selected.meta.txt",header=T)
microbio_selected <- microbio_selected[which(sapply(strsplit(as.character(microbio_selected$ID),"_"), `[`, 3)!="H2"),]
rownames(microbio_selected) <- sapply(strsplit(as.character(microbio_selected$ID),"_H"), `[`, 1)
microbio_selected$bmi_class <- factor(sapply(strsplit(as.character(microbio_selected$bmi_class),"-"), `[`, 1),levels = c("Lean","Overweight","Obese"))


######################
# FEATURES-COUNTS ----
######################

# Feature files: counts of library for each sample
# these files start with "#" which R does not read because interprets it as comment
# Read the header extract it and then name the columns from the file
ssu_feature_dada2 <- read.table("counts/16S_Dada2_Table.txt", header = F)
ssu_con <- file("counts/16S_Dada2_Table.txt", "r")
ssu_header <- readLines(ssu_con, n=2)
colnames(ssu_feature_dada2) <- strsplit(ssu_header, "\\s+")[[2]][-1]
close(ssu_con)

dnaK_feature_dada2 <- read.table("counts/dnaK_Dada2_Table.txt", header = F)
dnaK_con <- file("counts/dnaK_Dada2_Table.txt", "r")
dnaK_header_dada2 <- readLines(dnaK_con, n=2)
colnames(dnaK_feature_dada2) <- strsplit(dnaK_header_dada2, "\\s+")[[2]][-1]
close(dnaK_con)

gyrB_feature_dada2 <- read.table("counts/gyrB_Dada2_Table.txt", header = F)
gyrB_con <- file("counts/gyrB_Dada2_Table.txt", "r")
gyrB_header_dada2 <- readLines(gyrB_con, n=2)
colnames(gyrB_feature_dada2) <- strsplit(gyrB_header_dada2, "\\s+")[[2]][-1]   
close(gyrB_con)

# ASV table of dada2
ssu_mat_dada2 = as.matrix(ssu_feature_dada2[,-1])
rownames(ssu_mat_dada2) <- ssu_feature_dada2[,1]

dnaK_mat_dada2 = as.matrix(dnaK_feature_dada2[,-1])
rownames(dnaK_mat_dada2) <- dnaK_feature_dada2[,1]

gyrB_mat_dada2 = as.matrix(gyrB_feature_dada2[,-1])
rownames(gyrB_mat_dada2) <- gyrB_feature_dada2[,1]
                                                                                                                         
# OTU tables dada2
ssu_OTU_dada2 = otu_table(ssu_mat_dada2, taxa_are_rows = TRUE)
colnames(ssu_OTU_dada2) <- ssu_meta[match(colnames(ssu_OTU_dada2),ssu_meta$lib),]$V1
colnames(ssu_OTU_dada2) <- gsub("-", "_", colnames(ssu_OTU_dada2)) #change lib names for sample names (low stripes)
colnames(ssu_OTU_dada2) <- sapply(strsplit(as.character(colnames(ssu_OTU_dada2)),"_H"), `[`, 1)

dnaK_OTU_dada2 = otu_table(dnaK_mat_dada2, taxa_are_rows = TRUE)
colnames(dnaK_OTU_dada2) <- dnaK_meta[match(colnames(dnaK_OTU_dada2),dnaK_meta$lib),1] #change lib names for sample names
dnaK_OTU_dada2<-dnaK_OTU_dada2[,-which(colnames(dnaK_OTU_dada2)=="Bacterial_mix")] #remove bacterial mix

gyrB_OTU_dada2 = otu_table(gyrB_mat_dada2, taxa_are_rows = TRUE)
colnames(gyrB_OTU_dada2) <- gyrB_meta[match(colnames(gyrB_OTU_dada2),gyrB_meta$lib),1] #change lib names for sample names
gyrB_OTU_dada2<-gyrB_OTU_dada2[,-which(colnames(gyrB_OTU_dada2)=="Mix_bacteriano")] #remove bacterial mix
colnames(gyrB_OTU_dada2) <- sapply(strsplit(as.character(colnames(gyrB_OTU_dada2)),"_H"), `[`, 1)

###############
# TAXONOMY ----
###############

ssu_dada2_taxcomplete <- read.delim("taxonomy/16S_dada2_gtdb_taxonomy.tsv", sep="\t", header=T)
dnaK_dada2_taxcomplete <- read.csv("taxonomy/dnaK_dada2_gtdb_taxonomy.txt", header=T)
gyrB_dada2_taxcomplete <- read.csv("taxonomy/gyrB_dada2_gtdb_taxonomy.txt", header=T)

# Convert to matrix
ssu_taxmat_dada2 = as.matrix(ssu_dada2_taxcomplete)
rownames(ssu_taxmat_dada2) <- ssu_dada2_taxcomplete[,1]
ssu_TAX_dada2 = tax_table(ssu_taxmat_dada2) #[,-ncol(ssu_taxmat_dada2)]
colnames(ssu_TAX_dada2) <- c("seq_id","query","domain","phylum","class","order","family","genus","species")

dnaK_taxmat_dada2 = as.matrix(dnaK_dada2_taxcomplete)
rownames(dnaK_taxmat_dada2) <- dnaK_dada2_taxcomplete[,1]
dnaK_TAX_dada2 = tax_table(dnaK_taxmat_dada2) #[,-ncol(dnaK_taxmat_dada2)]
                                                                                                                         
gyrB_taxmat_dada2 = as.matrix(gyrB_dada2_taxcomplete)
rownames(gyrB_taxmat_dada2) <- gyrB_dada2_taxcomplete[,1]
gyrB_TAX_dada2 = tax_table(gyrB_taxmat_dada2)  #remove the last column of db -> NO necessary anymore!


############
# TREES ----
############

ssu_nwk_dada2 <- read.tree("trees/16S_rooted_tree_Dada2.nwk")
dnaK_qza_dada2 <- read_qza("trees/dnaK_rooted_tree_Dada2.qza")
gyrB_qza_dada2 <- read_qza("trees/gyrB_rooted_tree_Dada2.qza")


####################
# DNA SEQUENCES ----
####################

ssu_seq <- readDNAStringSet("sequences/16S_RepSeq_Dada2.fasta")
dnaK_seq <- readDNAStringSet("sequences/dnaK_Dada2_RepSeq.fasta")
gyrB_seq <- readDNAStringSet("sequences/gyrB_Dada2_RepSeq.fasta")


###############
# PHYLOSEQ ----
###############

ssu_physeq_dada2 = phyloseq(ssu_OTU_dada2, ssu_TAX_dada2, sample_data(microbio_selected), ssu_nwk_dada2, ssu_seq)
dnaK_physeq_dada2 = phyloseq(dnaK_OTU_dada2, dnaK_TAX_dada2, sample_data(microbio_selected), dnaK_qza_dada2$data, dnaK_seq)
gyrB_physeq_dada2 = phyloseq(gyrB_OTU_dada2, gyrB_TAX_dada2, sample_data(microbio_selected), gyrB_qza_dada2$data, gyrB_seq)

# Filter out ASVs that are only found once
ssu_physeq <- phyloseq::subset_samples(ssu_physeq_dada2, phyloseq::sample_sums(ssu_physeq_dada2) > 100)  
ssu_physeq <- phyloseq::prune_taxa(phyloseq::taxa_sums(ssu_physeq) > 2, ssu_physeq)

dnaK_physeq <- phyloseq::subset_samples(dnaK_physeq_dada2, phyloseq::sample_sums(dnaK_physeq_dada2) > 100)
dnaK_physeq <- phyloseq::prune_taxa(phyloseq::taxa_sums(dnaK_physeq) > 2, dnaK_physeq)

gyrB_physeq <- phyloseq::subset_samples(gyrB_physeq_dada2, phyloseq::sample_sums(gyrB_physeq_dada2) > 100)
gyrB_physeq <- phyloseq::prune_taxa(phyloseq::taxa_sums(gyrB_physeq) > 2, gyrB_physeq)

# Rarefaction
ssu_rare <- phyloseq::rarefy_even_depth(ssu_physeq, rngseed = 123, replace = FALSE)
dnaK_rare <- phyloseq::rarefy_even_depth(dnaK_physeq, rngseed = 123, replace = FALSE)
gyrB_rare <- phyloseq::rarefy_even_depth(gyrB_physeq, rngseed = 123, replace = FALSE)

                                                                                                                         
#############
# DESeq2 ----
#############

ssu_physeq <- prune_taxa(taxa_sums(ssu_physeq_dada2) > 1, ssu_physeq_dada2)
ssu_physeq <- prune_samples(sample_sums(ssu_physeq_dada2)>0, ssu_physeq_dada2)
ssu_dds = phyloseq_to_deseq2(ssu_physeq, ~ bmi_class)
ssu_dds.bmi <- estimateSizeFactors(ssu_dds,type="poscounts") # solves zero-inflated counts
ssu_dds.wald= DESeq(ssu_dds.bmi, fitType="local")

dnaK_physeq <- prune_taxa(taxa_sums(dnaK_physeq_dada2) > 1, dnaK_physeq_dada2)
dnaK_physeq <- prune_samples(sample_sums(dnaK_physeq_dada2)>0,dnaK_physeq_dada2)
dnaK_dds = phyloseq_to_deseq2(dnaK_physeq, ~ bmi_class)
dnaK_dds.bmi <- estimateSizeFactors(dnaK_dds,type="poscounts") # solves zero-inflated counts
dnaK_dds.wald= DESeq(dnaK_dds.bmi, fitType="local")

gyrB_physeq <- prune_taxa(taxa_sums(gyrB_physeq_dada2) > 1, gyrB_physeq_dada2)
gyrB_physeq <- prune_samples(sample_sums(gyrB_physeq_dada2)>0,gyrB_physeq_dada2)
gyrB_dds = phyloseq_to_deseq2(gyrB_physeq, ~ bmi_class)
gyrB_dds.bmi <- estimateSizeFactors(gyrB_dds,type="poscounts") # solves zero-inflated counts
gyrB_dds.wald <- DESeq(gyrB_dds.bmi, fitType="local")

                                                                                                                         
# Lean vs. obese
alpha = 0.05

ssu_res.nw_ob <- results(ssu_dds.wald,contrast=c("bmi_class","Obese","Lean"))
ssu_res.nw_ob = ssu_res.nw_ob[order(ssu_res.nw_ob$padj, na.last=NA), ]
ssu_sigtab.nw_ob = ssu_res.nw_ob[(ssu_res.nw_ob$padj < alpha), ]
ssu_sigtab.nw_ob = cbind(as(ssu_sigtab.nw_ob, "data.frame"), as(tax_table(ssu_physeq)[rownames(ssu_sigtab.nw_ob), ], "matrix"))

dnaK_res.nw_ob <- results(dnaK_dds.wald,contrast=c("bmi_class","Obese","Lean"))
dnaK_res.nw_ob = dnaK_res.nw_ob[order(dnaK_res.nw_ob$padj, na.last=NA), ]
dnaK_sigtab.nw_ob = dnaK_res.nw_ob[(dnaK_res.nw_ob$padj < alpha), ]
dnaK_sigtab.nw_ob = cbind(as(dnaK_sigtab.nw_ob, "data.frame"), as(tax_table(dnaK_physeq)[rownames(dnaK_sigtab.nw_ob), ], "matrix"))

gyrB_res.nw_ob <- results(gyrB_dds.wald,contrast=c("bmi_class","Obese","Lean"))
gyrB_res.nw_ob = gyrB_res.nw_ob[order(gyrB_res.nw_ob$padj, na.last=NA), ]
gyrB_sigtab.nw_ob = gyrB_res.nw_ob[(gyrB_res.nw_ob$padj < alpha), ]
gyrB_sigtab.nw_ob = cbind(as(gyrB_sigtab.nw_ob, "data.frame"), as(tax_table(gyrB_physeq)[rownames(gyrB_sigtab.nw_ob), ], "matrix"))
