#source("../scripts/phylotags_preprocessing.R")

# Test results when clustering at 97%
# 16S rRNA
ssu_refseq <- refseq(ssu_physeq_dada2) 
nproc <- 1 # Increase to use multiple processors
ssu_aln <- DECIPHER::AlignSeqs(ssu_refseq, processors = nproc)
ssu_d <- DECIPHER::DistanceMatrix(ssu_aln, processors = nproc)

# dnaK
dnaK_refseq <- refseq(dnaK_physeq_dada2)
nproc <- 1 # Increase to use multiple processors
dnaK_aln <- DECIPHER::AlignSeqs(dnaK_refseq, processors = nproc)
dnaK_d <- DECIPHER::DistanceMatrix(dnaK_aln, processors = nproc)

# gyrB 
gyrB_refseq <- refseq(gyrB_physeq_dada2)
nproc <- 1 # Increase to use multiple processors
gyrB_aln <- DECIPHER::AlignSeqs(gyrB_refseq, processors = nproc)
gyrB_d <- DECIPHER::DistanceMatrix(gyrB_aln, processors = nproc)


### Clusters
# 16S rRNA
ssu_clusters_97 <- DECIPHER::TreeLine(
  myDistMatrix=ssu_d, 
  method = "complete",
  cutoff = 0.03, # 0.1 corresponds to 90% OTUs, 0.16 to 84% (result from Figure 1)
  type = "clusters",
  processors = nproc
)

ssu_clusters_90 <- DECIPHER::TreeLine(
  myDistMatrix=ssu_d, 
  method = "complete",
  cutoff = 0.1, # 0.1 corresponds to 90% OTUs, 0.16 to 84% (result from Figure 1)
  type = "clusters",
  processors = nproc
)

ssu_clusters_85 <- DECIPHER::TreeLine(
  myDistMatrix=ssu_d, 
  method = "complete",
  cutoff = 0.15, # 0.1 corresponds to 90% OTUs, 0.16 to 84% (result from Figure 1)
  type = "clusters",
  processors = nproc
)


# dnaK
dnaK_clusters_97 <- DECIPHER::TreeLine(
  myDistMatrix=dnaK_d, 
  method = "complete",
  cutoff = 0.03, # 0.1 corresponds to 90% OTUs, 0.16 to 84% (result from Figure 1)
  type = "clusters",
  processors = nproc
)

dnaK_clusters_90 <- DECIPHER::TreeLine(
  myDistMatrix=dnaK_d, 
  method = "complete",
  cutoff = 0.1, # 0.1 corresponds to 90% OTUs, 0.16 to 84% (result from Figure 1)
  type = "clusters",
  processors = nproc
)

dnaK_clusters_85 <- DECIPHER::TreeLine(
  myDistMatrix=dnaK_d, 
  method = "complete",
  cutoff = 0.15, # 0.1 corresponds to 90% OTUs, 0.16 to 84% (result from Figure 1)
  type = "clusters",
  processors = nproc
)

# gyrB
gyrB_clusters_97 <- DECIPHER::TreeLine(
  myDistMatrix=gyrB_d, 
  method = "complete",
  cutoff = 0.03, # 0.1 corresponds to 90% OTUs, 0.16 to 84% (result from Figure 1)
  type = "clusters",
  processors = nproc
)

gyrB_clusters_90 <- DECIPHER::TreeLine(
  myDistMatrix=gyrB_d, 
  method = "complete",
  cutoff = 0.1, # 0.1 corresponds to 90% OTUs, 0.16 to 84% (result from Figure 1)
  type = "clusters",
  processors = nproc
)

gyrB_clusters_85 <- DECIPHER::TreeLine(
  myDistMatrix=gyrB_d, 
  method = "complete",
  cutoff = 0.15, # 0.1 corresponds to 90% OTUs, 0.16 to 84% (result from Figure 1)
  type = "clusters",
  processors = nproc
)

# Use speedyseq to merge taxa in phyloseq object, as explained in:
# https://github.com/mikemc/speedyseq/blob/main/NEWS.md#new-general-purpose-vectorized-merging-function 
library(speedyseq)

# 16S rRNA
ssu_physeq_clust85 <- merge_taxa_vec(ssu_physeq_dada2, group=ssu_clusters_85$cluster, tax_adjust = 0)
ssu_clusters_fulltax_clust85 <- merge(ssu_clusters_85,as.data.frame(ssu_TAX_dada2), by="row.names")
ssu_rare_85 <- phyloseq::rarefy_even_depth(ssu_physeq_clust85, rngseed = 123, replace = FALSE)

ssu_physeq_clust90 <- merge_taxa_vec(ssu_physeq_dada2, group=ssu_clusters_90$cluster, tax_adjust = 0)
ssu_clusters_fulltax_clust90 <- merge(ssu_clusters_90,as.data.frame(ssu_TAX_dada2), by="row.names")
ssu_rare_90 <- phyloseq::rarefy_even_depth(ssu_physeq_clust90, rngseed = 123, replace = FALSE)

ssu_physeq_clust97 <- merge_taxa_vec(ssu_physeq_dada2, group=ssu_clusters_97$cluster, tax_adjust = 0)
ssu_clusters_fulltax_clust97 <- merge(ssu_clusters_97,as.data.frame(ssu_TAX_dada2), by="row.names")
ssu_rare_97 <- phyloseq::rarefy_even_depth(ssu_physeq_clust97, rngseed = 123, replace = FALSE)

ssu_shannon.fam <- marker_tax_stats("ssu",85,"family")$shannon
ssu_shannon.gen <- marker_tax_stats("ssu",85,"genus")$shannon


# dnaK
dnaK_physeq_clust85 <- merge_taxa_vec(dnaK_physeq_dada2, group=dnaK_clusters_85$cluster, tax_adjust = 0)
dnaK_clusters_fulltax_clust85 <- merge(dnaK_clusters_85,as.data.frame(dnaK_TAX_dada2), by="row.names")
dnaK_rare_85 <- phyloseq::rarefy_even_depth(dnaK_physeq_clust85, rngseed = 123, replace = FALSE)

dnaK_physeq_clust90 <- merge_taxa_vec(dnaK_physeq_dada2, group=dnaK_clusters_90$cluster, tax_adjust = 0)
dnaK_clusters_fulltax_clust90 <- merge(dnaK_clusters_90,as.data.frame(dnaK_TAX_dada2), by="row.names")
dnaK_rare_90 <- phyloseq::rarefy_even_depth(dnaK_physeq_clust90, rngseed = 123, replace = FALSE)

dnaK_physeq_clust97 <- merge_taxa_vec(dnaK_physeq_dada2, group=dnaK_clusters_97$cluster, tax_adjust = 0)
dnaK_clusters_fulltax_clust97 <- merge(dnaK_clusters_97,as.data.frame(dnaK_TAX_dada2), by="row.names")
dnaK_rare_97 <- phyloseq::rarefy_even_depth(dnaK_physeq_clust97, rngseed = 123, replace = FALSE)

# gyrB
gyrB_physeq_clust85 <- merge_taxa_vec(gyrB_physeq_dada2, group=gyrB_clusters_85$cluster, tax_adjust = 0)
gyrB_clusters_fulltax_clust85 <- merge(gyrB_clusters_85,as.data.frame(gyrB_TAX_dada2), by="row.names")
gyrB_rare_85 <- phyloseq::rarefy_even_depth(gyrB_physeq_clust85, rngseed = 123, replace = FALSE, sample.size = 82)

gyrB_physeq_clust90 <- merge_taxa_vec(gyrB_physeq_dada2, group=gyrB_clusters_90$cluster, tax_adjust = 0)
gyrB_clusters_fulltax_clust90 <- merge(gyrB_clusters_90,as.data.frame(gyrB_TAX_dada2), by="row.names")
gyrB_rare_90 <- phyloseq::rarefy_even_depth(gyrB_physeq_clust90, rngseed = 123, replace = FALSE, sample.size = 82)

gyrB_physeq_clust97 <- merge_taxa_vec(gyrB_physeq_dada2, group=gyrB_clusters_97$cluster, tax_adjust = 0)
gyrB_clusters_fulltax_clust97 <- merge(gyrB_clusters_97,as.data.frame(gyrB_TAX_dada2), by="row.names")
#pruned <- prune_samples(sample_sums(gyrB_clusters_fulltax_clust97) > 0, gyrB_clusters_fulltax_clust97,taxa_are_rows=TRUE)
gyrB_rare_97 <- phyloseq::rarefy_even_depth(gyrB_physeq_clust97, rngseed = 123, replace = FALSE, sample.size = 82)

gyrB_shannon.fam <- marker_tax_stats("gyrB",97,"family")$shannon
gyrB_shannon.gen <- marker_tax_stats("gyrB",97,"genus")$shannon

gyrB_shannon.fam <- marker_tax_stats("gyrB",85,"family")$shannon
gyrB_shannon.gen <- marker_tax_stats("gyrB",85,"genus")$shannon


#########################################
# Group at different identity levels ----
#########################################

# Sequence identity = 97%
rm(ssu_rare, dnaK_rare, gyrB_rare)
idlevel=97
ssu_rare <- ssu_rare_97
dnaK_rare <- dnaK_rare_97
gyrB_rare <- gyrB_rare_97
# From here on, continue with the code of phylotags_diversity.R script

# Sequence identity = 90%
rm(ssu_rare, dnaK_rare, gyrB_rare)
idlevel=90
ssu_rare <- ssu_rare_90
dnaK_rare <- dnaK_rare_90
gyrB_rare <- gyrB_rare_90
# From here on, continue with the code of phylotags_diversity.R script

# Sequence identity = 85%
rm(ssu_rare, dnaK_rare, gyrB_rare)
idlevel=85
ssu_rare <- ssu_rare_85
dnaK_rare <- dnaK_rare_85
gyrB_rare <- gyrB_rare_85
# From here on, continue with the code of phylotags_diversity.R script
