##################################
# Select microbes of interest ----
##################################

# From a UNIX prompt, do this (modify taxa accordingly):
#grep -w "g__CAG-83" Global_metagenomes_counts.tsv > Global_metagenomes_counts.CAG_83.tsv

###############
# METADATA ----
###############

global_metadata <- read.table("global_metagenomes/Global_metagenomes_metadata.tsv", header = T, sep = "\t")

# Extract samples with BMI annotation (drop NAs)
global_metadata <- global_metadata %>% drop_na("BMI")

# Select samples with disease annotations of interest; keep NAs here (replaced by unkown)
global_metadata <- global_metadata %>% 
  replace_na(list(disease = "unknown")) %>%
  filter(disease!="fatty_liver;hypertension" & disease!="fatty_liver" & disease!="fatty_liver;T2D;hypertension" & disease!="fatty_liver;T2D")

# Define the row names
global_metadata <- column_to_rownames(global_metadata, var = "Sample")
global_metadata <- global_metadata[order(rownames(global_metadata)),]

# Remove individuals with low BMI (<18.5 kg/m2)
global_metadata <- global_metadata[global_metadata$BMI>=18.5,]

# BMI categories
global_metadata$bmi_class <- ifelse(global_metadata$BMI>=18.5 & global_metadata$BMI<25, "Lean",
                                    ifelse(global_metadata$BMI>=25 & global_metadata$BMI<30, "Overweight",
                                                  "Obese"))


######################
# FEATURES-COUNTS ----
######################

# Since files are large, gunzip them befre runnig the following chunk
# From a UNIX prompt, do this:
#gunzip *.gz

# Oscillospiraceae
oscillo_counts <- read.table("global_metagenomes/Global_metagenomes_counts.Ocillospiraceae.tsv", header = T, sep = "\t")
oscillo_sample <- unique(oscillo_counts$Sample)

# Trim samples to those previously selected in global_metadata
oscillo_counts <- oscillo_counts[oscillo_sample %in% rownames(global_metadata),]

# Create an OTU table ("unmelt" the counts dataset)
oscillo_otu <- oscillo_counts[order(oscillo_counts$Sample, oscillo_counts$name),]
oscillo_otu <- reshape2::dcast(data = oscillo_otu,formula = Sample~name, fun.aggregate = sum, value.var = "Abundance")

# Define the row names from the otu column and remove the column since it is now used as a row name
oscillo_otu <- column_to_rownames(oscillo_otu, var = "Sample")

# Verify that the samples were correctly selected
identical(rownames(oscillo_otu), rownames(global_metadata))


# Acutalibacteraceae
acutali_counts <- read.table("global_metagenomes/Global_metagenomes_counts.Acutalibacteraceae.tsv", header = T, sep = "\t")
acutali_sample <- unique(acutali_counts$Sample)

# Trim samples to those previously selected in global_metadata
acutali_counts <- acutali_counts[acutali_sample %in% rownames(global_metadata),]

# Create an OTU table ("unmelt" the counts dataset)
acutali_otu <- acutali_counts[order(acutali_counts$Sample, acutali_counts$name),]
acutali_otu <- reshape2::dcast(data = acutali_otu,formula = Sample~name, fun.aggregate = sum, value.var = "Abundance")

# Define the row names from the otu column and remove the column since it is now used as a row name
acutali_otu <- column_to_rownames(acutali_otu, var = "Sample")

# Verify that the samples were correctly selected
identical(rownames(acutali_otu), rownames(global_metadata))


# Ruminococcaceae
rumino_counts <- read.table("global_metagenomes/Global_metagenomes_counts.Ruminococcaceae.tsv", header = T, sep = "\t")
rumino_sample <- unique(rumino_counts$Sample)

# Trim samples to those previously selected in global_metadata
rumino_counts <- rumino_counts[rumino_sample %in% rownames(global_metadata),]

# Create an OTU table ("unmelt" the counts dataset)
rumino_otu <- rumino_counts[order(rumino_counts$Sample, rumino_counts$name),]
rumino_otu <- reshape2::dcast(data = rumino_otu,formula = Sample~name, fun.aggregate = sum, value.var = "Abundance")

# Define the row names from the otu column and remove the column since it is now used as a row name
rumino_otu <- column_to_rownames(rumino_otu, var = "Sample")

# Verify that the samples were correctly selected
identical(rownames(rumino_otu), rownames(global_metadata))


# Lachnospiraceae
lachno_counts <- read.table("global_metagenomes/Global_metagenomes_counts.Lachnospiraceae.tsv", header = T, sep = "\t")
lachno_sample <- unique(lachno_counts$Sample)

# Trim samples to those previously selected in global_metadata
lachno_counts <- lachno_counts[lachno_sample %in% rownames(global_metadata),]

# Create an OTU table ("unmelt" the counts dataset)
lachno_otu <- lachno_counts[order(lachno_counts$Sample, lachno_counts$name),]
lachno_otu <- reshape2::dcast(data = lachno_otu,formula = Sample~name, fun.aggregate = sum, value.var = "Abundance")

# Define the row names from the otu column and remove the column since it is now used as a row name
lachno_otu <- column_to_rownames(lachno_otu, var = "Sample")

# Verify that the samples were correctly selected
identical(rownames(lachno_otu), rownames(global_metadata))


# g__CAG-83 (Oscillospiraceae)
cag83_counts <- read.table("global_metagenomes/Global_metagenomes_counts.CAG_83.tsv", header = T, sep = "\t")
cag83_sample <- unique(cag83_counts$Sample)

# Trim samples to those previously selected in global_metadata
cag83_counts <- cag83_counts[cag83_sample %in% rownames(global_metadata),]

# Create an OTU table ("unmelt" the counts dataset)
cag83_otu <- cag83_counts[order(cag83_counts$Sample, cag83_counts$name),]
cag83_otu <- reshape2::dcast(data = cag83_otu,formula = Sample~name, fun.aggregate = sum, value.var = "Abundance")

# Define the row names from the otu column and remove the column since it is now used as a row name
cag83_otu <- column_to_rownames(cag83_otu, var = "Sample")

# Verify that the samples were correctly selected
identical(rownames(cag83_otu), rownames(global_metadata))


# g__CAG-170 (Oscillospiraceae)
cag170_counts <- read.table("global_metagenomes/Global_metagenomes_counts.CAG_170.tsv", header = T, sep = "\t")
cag170_sample <- unique(cag170_counts$Sample)

# Trim samples to those previously selected in global_metadata
cag170_counts <- cag170_counts[cag170_sample %in% rownames(global_metadata),]

# Create an OTU table ("unmelt" the counts dataset)
cag170_otu <- cag170_counts[order(cag170_counts$Sample, cag170_counts$name),]
cag170_otu <- reshape2::dcast(data = cag170_otu,formula = Sample~name, fun.aggregate = sum, value.var = "Abundance")

# Define the row names from the otu column and remove the column since it is now used as a row name
cag170_otu <- column_to_rownames(cag170_otu, var = "Sample")

# Verify that the samples were correctly selected
identical(rownames(cag170_otu), rownames(global_metadata))


# g__CAG-177 (Acutalibacteraceae)
cag177_counts <- read.table("global_metagenomes/Global_metagenomes_counts.CAG_177.tsv", header = T, sep = "\t")
cag177_sample <- unique(cag177_counts$Sample)

# Trim samples to those previously selected in global_metadata
cag177_counts <- cag177_counts[cag177_sample %in% rownames(global_metadata),]

# Create an OTU table ("unmelt" the counts dataset)
cag177_otu <- cag177_counts[order(cag177_counts$Sample, cag177_counts$name),]
cag177_otu <- reshape2::dcast(data = cag177_otu,formula = Sample~name, fun.aggregate = sum, value.var = "Abundance")

# Define the row names from the otu column and remove the column since it is now used as a row name
cag177_otu <- column_to_rownames(cag177_otu, var = "Sample")

# Verify that the samples were correctly selected
identical(rownames(cag177_otu), rownames(global_metadata))


# g__CAG-81 (Lachnospiraceae)
cag81_counts <- read.table("global_metagenomes/Global_metagenomes_counts.CAG_81.tsv", header = T, sep = "\t")
cag81_sample <- unique(cag81_counts$Sample)

# Trim samples to those previously selected in global_metadata
cag81_counts <- cag81_counts[cag81_sample %in% rownames(global_metadata),]

# Create an OTU table ("unmelt" the counts dataset)
cag81_otu <- cag81_counts[order(cag81_counts$Sample, cag81_counts$name),]
cag81_otu <- reshape2::dcast(data = cag81_otu,formula = Sample~name, fun.aggregate = sum, value.var = "Abundance")

# Define the row names from the otu column and remove the column since it is now used as a row name
cag81_otu <- column_to_rownames(cag81_otu, var = "Sample")

# Verify that the samples were correctly selected
identical(rownames(cag81_otu), rownames(global_metadata))


### g__COE1 (Lachnospiraceae)
coe1_counts <- read.table("global_metagenomes/Global_metagenomes_counts.COE1.tsv", header = T, sep = "\t")
coe1_sample <- unique(coe1_counts$Sample)

# Trim samples to those previously selected in global_metadata
coe1_counts <- coe1_counts[coe1_sample %in% rownames(global_metadata),]

# Create an OTU table ("unmelt" the counts dataset)
coe1_otu <- coe1_counts[order(coe1_counts$Sample, coe1_counts$name),]
coe1_otu <- reshape2::dcast(data = coe1_otu,formula = Sample~name, fun.aggregate = sum, value.var = "Abundance")

# Define the row names from the otu column and remove the column since it is now used as a row name
coe1_otu <- column_to_rownames(coe1_otu, var = "Sample")

# Verify that the samples were correctly selected
identical(rownames(coe1_otu), rownames(global_metadata))


# g__Ruminococcus_A (Lachnospiraceae)
ruminoA_counts <- read.table("global_metagenomes/Global_metagenomes_counts.Ruminococcus_A.tsv", header = T, sep = "\t")
ruminoA_sample <- unique(ruminoA_counts$Sample)

# Trim samples to those previously selected in global_metadata
ruminoA_counts <- ruminoA_counts[ruminoA_sample %in% rownames(global_metadata),]

# Create an OTU table ("unmelt" the counts dataset)
ruminoA_otu <- ruminoA_counts[order(ruminoA_counts$Sample, ruminoA_counts$name),]
ruminoA_otu <- reshape2::dcast(data = ruminoA_otu,formula = Sample~name, fun.aggregate = sum, value.var = "Abundance")

# Define the row names from the otu column and remove the column since it is now used as a row name
ruminoA_otu <- column_to_rownames(ruminoA_otu, var = "Sample")

# Verify that the samples were correctly selected
identical(rownames(ruminoA_otu), rownames(global_metadata))


###############
# TAXONOMY ----
###############

# Oscillospiraceae
oscillo_tax <- data.frame(taxonomy=oscillo_counts$taxonomy)
oscillo_tax <- unique(oscillo_tax)
oscillo_tax <- oscillo_tax %>% 
  separate(taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")
rownames(oscillo_tax) <- oscillo_tax$Species

# Verify that the species were correctly selected
identical(sort(colnames(oscillo_otu)), sort(rownames(oscillo_tax)))


# Ruminococcaceae
rumino_tax <- data.frame(taxonomy=rumino_counts$taxonomy)
rumino_tax <- unique(rumino_tax)
rumino_tax <- rumino_tax %>% 
  separate(taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")
rownames(rumino_tax) <- rumino_tax$Species

# Verify that the species were correctly selected
identical(sort(colnames(rumino_otu)), sort(rownames(rumino_tax)))


# Acutalibacteraceae
acutali_tax <- data.frame(taxonomy=acutali_counts$taxonomy)
acutali_tax <- unique(acutali_tax)
acutali_tax <- acutali_tax %>% 
  separate(taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")
rownames(acutali_tax) <- acutali_tax$Species

# Verify that the species were correctly selected
identical(sort(colnames(acutali_otu)), sort(rownames(acutali_tax)))


# Ruminococcaceae
rumino_tax <- data.frame(taxonomy=rumino_counts$taxonomy)
rumino_tax <- unique(rumino_tax)
rumino_tax <- rumino_tax %>% 
  separate(taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")
rownames(rumino_tax) <- rumino_tax$Species

# Verify that the species were correctly selected
identical(sort(colnames(rumino_otu)), sort(rownames(rumino_tax)))


# Lachnospiraceae
lachno_tax <- data.frame(taxonomy=lachno_counts$taxonomy)
lachno_tax <- unique(lachno_tax)
lachno_tax <- lachno_tax %>% 
  separate(taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")
rownames(lachno_tax) <- lachno_tax$Species

# Verify that the species were correctly selected
identical(sort(colnames(lachno_otu)), sort(rownames(lachno_tax)))


# g__CAG-83 (Oscillospiraceae)
cag83_tax <- data.frame(taxonomy=cag83_counts$taxonomy)
cag83_tax <- unique(cag83_tax)
cag83_tax <- cag83_tax %>% 
  separate(taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")
rownames(cag83_tax) <- cag83_tax$Species

# Verify that the species were correctly selected
identical(sort(colnames(cag83_otu)), sort(rownames(cag83_tax)))


# g__ CAG-170 (Oscillospiraceae)
cag170_tax <- data.frame(taxonomy=cag170_counts$taxonomy)
cag170_tax <- unique(cag170_tax)
cag170_tax <- cag170_tax %>% 
  separate(taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")
rownames(cag170_tax) <- cag170_tax$Species

# Verify that the species were correctly selected
identical(sort(colnames(cag170_otu)), sort(rownames(cag170_tax)))


# g__CAG-177 (Acutalibacteraceae)
cag177_tax <- data.frame(taxonomy=cag177_counts$taxonomy)
cag177_tax <- unique(cag177_tax)
cag177_tax <- cag177_tax %>% 
  separate(taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")
rownames(cag177_tax) <- cag177_tax$Species

# Verify that the species were correctly selected
identical(sort(colnames(cag177_otu)), sort(rownames(cag177_tax)))


# g__CAG-81 (Lachnospiraceae)
cag81_tax <- data.frame(taxonomy=cag81_counts$taxonomy)
cag81_tax <- unique(cag81_tax)
cag81_tax <- cag81_tax %>% 
  separate(taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")
rownames(cag81_tax) <- cag81_tax$Species

# Verify that the species were correctly selected
identical(sort(colnames(cag81_otu)), sort(rownames(cag81_tax)))


# g__COE1 (Lachnospiraceae)
coe1_tax <- data.frame(taxonomy=coe1_counts$taxonomy)
coe1_tax <- unique(coe1_tax)
coe1_tax <- coe1_tax %>% 
  separate(taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")
rownames(coe1_tax) <- coe1_tax$Species

# Verify that the species were correctly selected
identical(sort(colnames(coe1_otu)), sort(rownames(coe1_tax)))


# g__Ruminococcus_A (Lachnospiraceae)
ruminoA_tax <- data.frame(taxonomy=ruminoA_counts$taxonomy)
ruminoA_tax <- unique(ruminoA_tax)
ruminoA_tax <- ruminoA_tax %>% 
  separate(taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")
rownames(ruminoA_tax) <- ruminoA_tax$Species

# Verify that the species were correctly selected
identical(sort(colnames(ruminoA_otu)), sort(rownames(ruminoA_tax)))


###############
# PHYLOSEQ ----
###############

### Oscillospiraceae
# Transform into matrices OTU and taxonomy tables 
oscillo_otu <- as.matrix(oscillo_otu)
oscillo_tax <- as.matrix(oscillo_tax)

# phyloseq object
oscillo_physeq <- phyloseq(otu_table(oscillo_otu, taxa_are_rows = FALSE), tax_table(oscillo_tax), sample_data(global_metadata))


### Ruminococcaceae
# Transform into matrices OTU and taxonomy tables 
rumino_otu <- as.matrix(rumino_otu)
rumino_tax <- as.matrix(rumino_tax)

# phyloseq object
rumino_physeq <- phyloseq(otu_table(rumino_otu, taxa_are_rows = FALSE), tax_table(rumino_tax), sample_data(global_metadata))


### Acutalibacteraceae
# Transform into matrices OTU and taxonomy tables 
acutali_otu <- as.matrix(acutali_otu)
acutali_tax <- as.matrix(acutali_tax)

# phyloseq object
acutali_physeq <- phyloseq(otu_table(acutali_otu, taxa_are_rows = FALSE), tax_table(acutali_tax), sample_data(global_metadata))


### Lachnospiraceae
# Transform into matrices OTU and taxonomy tables 
lachno_otu <- as.matrix(lachno_otu)
lachno_tax <- as.matrix(lachno_tax)

# phyloseq object
lachno_physeq <- phyloseq(otu_table(lachno_otu, taxa_are_rows = FALSE), tax_table(lachno_tax), sample_data(global_metadata))


### g__CAG-83 (Oscillospiraceae)
# Transform into matrices OTU and taxonomy tables 
cag83_otu <- as.matrix(cag83_otu)
cag83_tax <- as.matrix(cag83_tax)

# phyloseq object
cag83_physeq <- phyloseq(otu_table(cag83_otu, taxa_are_rows = FALSE), tax_table(cag83_tax), sample_data(global_metadata))


### g__CAG-170 (Oscillospiraceae)
# Transform into matrices OTU and taxonomy tables 
cag170_otu <- as.matrix(cag170_otu)
cag170_tax <- as.matrix(cag170_tax)

# phyloseq object
cag170_physeq <- phyloseq(otu_table(cag170_otu, taxa_are_rows = FALSE), tax_table(cag170_tax), sample_data(global_metadata))


### g__CAG-177 (Acutalibacteraceae)
# Transform into matrices OTU and taxonomy tables 
cag177_otu <- as.matrix(cag177_otu)
cag177_tax <- as.matrix(cag177_tax)

# phyloseq object
cag177_physeq <- phyloseq(otu_table(cag177_otu, taxa_are_rows = FALSE), tax_table(cag177_tax), sample_data(global_metadata))


### g__CAG-81 (Lachnospiraceae)
# Transform into matrices OTU and taxonomy tables 
cag81_otu <- as.matrix(cag81_otu)
cag81_tax <- as.matrix(cag81_tax)

# phyloseq object
cag81_physeq <- phyloseq(otu_table(cag81_otu, taxa_are_rows = FALSE), tax_table(cag81_tax), sample_data(global_metadata))


### g__COE1 (Lachnospiraceae)
# Transform into matrices OTU and taxonomy tables 
coe1_otu <- as.matrix(coe1_otu)
coe1_tax <- as.matrix(coe1_tax)

# phyloseq object
coe1_physeq <- phyloseq(otu_table(coe1_otu, taxa_are_rows = FALSE), tax_table(coe1_tax), sample_data(global_metadata))


### g__Ruminococcus_A (Lachnospiraceae)
# Transform into matrices OTU and taxonomy tables 
ruminoA_otu <- as.matrix(ruminoA_otu)
ruminoA_tax <- as.matrix(ruminoA_tax)

# phyloseq object
ruminoA_physeq <- phyloseq(otu_table(ruminoA_otu, taxa_are_rows = FALSE), tax_table(ruminoA_tax), sample_data(global_metadata))


###########################################
# ALPHA DIVERSITY AND BMI ASSOCIATIONS ----
###########################################

### Oscillospiraceae
# Calculate alpha diversity for Oscillospiraceae
oscillo_global_adiv <- data.frame(
  "Observed" = phyloseq::estimate_richness(oscillo_physeq, measures = "Observed"),
  "shannon" = phyloseq::estimate_richness(oscillo_physeq, measures = "Shannon"),
  "shannon_bin" = ifelse(phyloseq::estimate_richness(oscillo_physeq, measures = "Shannon")!=0,1,0),
  "BMI" = phyloseq::sample_data(oscillo_physeq)$BMI,
  "bmi_class" = phyloseq::sample_data(oscillo_physeq)$bmi_class,
  "country" = phyloseq::sample_data(oscillo_physeq)$country,
  "dataset_name" = phyloseq::sample_data(oscillo_physeq)$dataset_name,
  "sex" = phyloseq::sample_data(oscillo_physeq)$gender,
  "age" = phyloseq::sample_data(oscillo_physeq)$age,
  "non_westernized" = phyloseq::sample_data(oscillo_physeq)$non_westernized)
colnames(oscillo_global_adiv) <- c("Observed", "shannon", "shannon_bin", "bmi", "bmi_class", "country", "dataset_name", "sex", "age", "non_westernized")

# BINARY MODEL (glm.bin)
Anova(glm(shannon_bin ~ bmi + sex + age, data = oscillo_global_adiv, family = binomial(link = logit)))

# CONTINUOUS MODEL (glm.cont)
# remove columns where all values are 0 (for numeric values)   
oscillo_shannon_cont <- oscillo_global_adiv[oscillo_global_adiv$shannon >0,]
Anova(glm(shannon ~ bmi + sex + age, data = oscillo_shannon_cont, family = Gamma(link = log)))

# scatterplot
ggplot(oscillo_shannon_cont, aes(x = bmi, y = shannon, color = bmi_class, shape = non_westernized)) +
  geom_point()


### Ruminococcaceae
# Calculate alpha diversity for Oscillospiraceae
rumino_global_adiv <- data.frame(
  "Observed" = phyloseq::estimate_richness(rumino_physeq, measures = "Observed"),
  "shannon" = phyloseq::estimate_richness(rumino_physeq, measures = "Shannon"),
  "shannon_bin" = ifelse(phyloseq::estimate_richness(rumino_physeq, measures = "Shannon")!=0,1,0),
  "BMI" = phyloseq::sample_data(rumino_physeq)$BMI,
  "bmi_class" = phyloseq::sample_data(rumino_physeq)$bmi_class,
  "country" = phyloseq::sample_data(rumino_physeq)$country,
  "dataset_name" = phyloseq::sample_data(rumino_physeq)$dataset_name,
  "sex" = phyloseq::sample_data(rumino_physeq)$gender,
  "age" = phyloseq::sample_data(rumino_physeq)$age,
  "non_westernized" = phyloseq::sample_data(rumino_physeq)$non_westernized)
colnames(rumino_global_adiv) <- c("Observed", "shannon", "shannon_bin", "bmi", "bmi_class", "country", "dataset_name", "sex", "age", "non_westernized")

# BINARY MODEL (glm.bin)
Anova(glm(shannon_bin ~ bmi + sex + age, data = rumino_global_adiv, family = binomial(link = logit)))

# CONTINUOUS MODEL (glm.cont)
# remove columns where all values are 0 (for numeric values)   
rumino_shannon_cont <- rumino_global_adiv[rumino_global_adiv$shannon >0,]
Anova(glm(shannon ~ bmi + sex + age, data = rumino_shannon_cont, family = Gamma(link = log)))

# scatterplot
ggplot(rumino_shannon_cont, aes(x = bmi, y = shannon, color = bmi_class, shape = non_westernized)) +
  geom_point()


### Acutalibacteraceae
# Calculate alpha diversity for Oscillospiraceae
acutali_global_adiv <- data.frame(
  "Observed" = phyloseq::estimate_richness(acutali_physeq, measures = "Observed"),
  "shannon" = phyloseq::estimate_richness(acutali_physeq, measures = "Shannon"),
  "shannon_bin" = ifelse(phyloseq::estimate_richness(acutali_physeq, measures = "Shannon")!=0,1,0),
  "BMI" = phyloseq::sample_data(acutali_physeq)$BMI,
  "bmi_class" = phyloseq::sample_data(acutali_physeq)$bmi_class,
  "country" = phyloseq::sample_data(acutali_physeq)$country,
  "dataset_name" = phyloseq::sample_data(acutali_physeq)$dataset_name,
  "sex" = phyloseq::sample_data(acutali_physeq)$gender,
  "age" = phyloseq::sample_data(acutali_physeq)$age,
  "non_westernized" = phyloseq::sample_data(acutali_physeq)$non_westernized)
colnames(acutali_global_adiv) <- c("Observed", "shannon", "shannon_bin", "bmi", "bmi_class", "country", "dataset_name", "sex", "age", "non_westernized")

# BINARY MODEL (glm.bin)
Anova(glm(shannon_bin ~ bmi + sex + age, data = acutali_global_adiv, family = binomial(link = logit)))

# CONTINUOUS MODEL
# remove columns where all values are 0 (for numeric values)   
acutali_shannon_cont <- acutali_global_adiv[acutali_global_adiv$shannon >0,]
Anova(glm(shannon ~ bmi + sex + age, data = acutali_shannon_cont, family = Gamma(link = log)))

# scatterplot
ggplot(acutali_shannon_cont, aes(x = bmi, y = shannon, color = bmi_class, shape = non_westernized)) +
  geom_point()


### Lachnospiraceae
# Calculate alpha diversity for Oscillospiraceae
lachno_global_adiv <- data.frame(
  "Observed" = phyloseq::estimate_richness(lachno_physeq, measures = "Observed"),
  "shannon" = phyloseq::estimate_richness(lachno_physeq, measures = "Shannon"),
  "shannon_bin" = ifelse(phyloseq::estimate_richness(lachno_physeq, measures = "Shannon")!=0,1,0),
  "BMI" = phyloseq::sample_data(lachno_physeq)$BMI,
  "bmi_class" = phyloseq::sample_data(lachno_physeq)$bmi_class,
  "country" = phyloseq::sample_data(lachno_physeq)$country,
  "dataset_name" = phyloseq::sample_data(lachno_physeq)$dataset_name,
  "sex" = phyloseq::sample_data(lachno_physeq)$gender,
  "age" = phyloseq::sample_data(lachno_physeq)$age,
  "non_westernized" = phyloseq::sample_data(lachno_physeq)$non_westernized)
colnames(lachno_global_adiv) <- c("Observed", "shannon", "shannon_bin", "bmi", "bmi_class", "country", "dataset_name", "sex", "age", "non_westernized")

# BINARY MODEL (glm.bin)
Anova(glm(shannon_bin ~ bmi + sex + age, data = lachno_global_adiv, family = binomial(link = logit)))

# CONTINUOUS MODEL (glm.cont)
# remove columns where all values are 0 (for numeric values)   
lachno_shannon_cont <- lachno_global_adiv[lachno_global_adiv$shannon >0,]
Anova(glm(shannon ~ bmi + sex + age, data = lachno_shannon_cont, family = Gamma(link = log)))

# scatterplot
ggplot(lachno_shannon_cont, aes(x = bmi, y = shannon, color = bmi_class, shape = non_westernized)) +
  geom_point()


### g__CAG-83 (Oscillospiraceae)
# Calculate alpha diversity for CAG-83
cag83_global_adiv <- data.frame(
  "Observed" = phyloseq::estimate_richness(cag83_physeq, measures = "Observed"),
  "shannon" = phyloseq::estimate_richness(cag83_physeq, measures = "Shannon"),
  "shannon_bin" = ifelse(phyloseq::estimate_richness(cag83_physeq, measures = "Shannon")!=0,1,0),
  "BMI" = phyloseq::sample_data(cag83_physeq)$BMI,
  "bmi_class" = phyloseq::sample_data(cag83_physeq)$bmi_class,
  "country" = phyloseq::sample_data(cag83_physeq)$country,
  "dataset_name" = phyloseq::sample_data(cag83_physeq)$dataset_name,
  "sex" = phyloseq::sample_data(cag83_physeq)$gender,
  "age" = phyloseq::sample_data(cag83_physeq)$age,
  "non_westernized" = phyloseq::sample_data(cag83_physeq)$non_westernized)
colnames(cag83_global_adiv) <- c("Observed", "shannon", "shannon_bin", "bmi", "bmi_class", "country", "dataset_name", "sex", "age", "non_westernized")

# BINARY MODEL (glm.bin)
Anova(glm(shannon_bin ~ bmi + sex + age, data = cag83_global_adiv, family = binomial(link = logit)))

# CONTINUOUS MODEL (glm.cont)
# remove columns where all values are 0 (for numeric values)   
cag83_shannon_cont <- cag83_global_adiv[cag83_global_adiv$shannon >0,]
Anova(glm(shannon ~ bmi + sex + age, data = cag83_shannon_cont, family = Gamma(link = log)))

# scatterplot
ggplot(cag83_shannon_cont, aes(x = bmi, y = shannon, color = bmi_class, shape = non_westernized)) +
  geom_point()


### g__CAG-170 (Oscillospiraceae)
# Calculate alpha diversity for CAG-170
cag170_global_adiv <- data.frame(
  "Observed" = phyloseq::estimate_richness(cag170_physeq, measures = "Observed"),
  "shannon" = phyloseq::estimate_richness(cag170_physeq, measures = "Shannon"),
  "shannon_bin" = ifelse(phyloseq::estimate_richness(cag170_physeq, measures = "Shannon")!=0,1,0),
  "BMI" = phyloseq::sample_data(cag170_physeq)$BMI,
  "bmi_class" = phyloseq::sample_data(cag170_physeq)$bmi_class,
  "country" = phyloseq::sample_data(cag170_physeq)$country,
  "dataset_name" = phyloseq::sample_data(cag170_physeq)$dataset_name,
  "sex" = phyloseq::sample_data(cag170_physeq)$gender,
  "age" = phyloseq::sample_data(cag170_physeq)$age,
  "non_westernized" = phyloseq::sample_data(cag170_physeq)$non_westernized)
colnames(cag170_global_adiv) <- c("Observed", "shannon", "shannon_bin", "bmi", "bmi_class", "country", "dataset_name", "sex", "age", "non_westernized")

# BINARY MODEL (glm.bin)
Anova(glm(shannon_bin ~ bmi + sex + age, data = cag170_global_adiv, family = binomial(link = logit)))

# CONTINUOUS MODEL (glm.cont)
# remove columns where all values are 0 (for numeric values)   
cag170_shannon_cont <- cag170_global_adiv[cag170_global_adiv$shannon >0,]
Anova(glm(shannon ~ bmi + sex + age, data = cag170_shannon_cont, family = Gamma(link = log)))

# scatterplot
ggplot(cag170_shannon_cont, aes(x = bmi, y = shannon, color = bmi_class, shape = non_westernized)) +
  geom_point()


### g__CAG-177 (Acutalibacteraceae)
# Calculate alpha diversity for CAG-177
cag177_global_adiv <- data.frame(
  "Observed" = phyloseq::estimate_richness(cag177_physeq, measures = "Observed"),
  "shannon" = phyloseq::estimate_richness(cag177_physeq, measures = "Shannon"),
  "shannon_bin" = ifelse(phyloseq::estimate_richness(cag177_physeq, measures = "Shannon")!=0,1,0),
  "BMI" = phyloseq::sample_data(cag177_physeq)$BMI,
  "bmi_class" = phyloseq::sample_data(cag177_physeq)$bmi_class,
  "country" = phyloseq::sample_data(cag177_physeq)$country,
  "dataset_name" = phyloseq::sample_data(cag177_physeq)$dataset_name,
  "sex" = phyloseq::sample_data(cag177_physeq)$gender,
  "age" = phyloseq::sample_data(cag177_physeq)$age,
  "non_westernized" = phyloseq::sample_data(cag177_physeq)$non_westernized)
colnames(cag177_global_adiv) <- c("Observed", "shannon", "shannon_bin", "bmi", "bmi_class", "country", "dataset_name", "sex", "age", "non_westernized")

# BINARY MODEL (glm.bin)
Anova(glm(shannon_bin ~ bmi + sex + age, data = cag177_global_adiv, family = binomial(link = logit)))

# CONTINUOUS MODEL (glm.cont)
# remove columns where all values are 0 (for numeric values)   
cag177_shannon_cont <- cag177_global_adiv[cag177_global_adiv$shannon >0,]
Anova(glm(shannon ~ bmi + sex + age, data = cag177_shannon_cont, family = Gamma(link = log)))

# scatterplot
ggplot(cag177_shannon_cont, aes(x = bmi, y = shannon, color = bmi_class, shape = non_westernized)) +
  geom_point()


### g__CAG-81 (Lachnospiraceae)
# Calculate alpha diversity for CAG-81
cag81_global_adiv <- data.frame(
  "Observed" = phyloseq::estimate_richness(cag81_physeq, measures = "Observed"),
  "shannon" = phyloseq::estimate_richness(cag81_physeq, measures = "Shannon"),
  "shannon_bin" = ifelse(phyloseq::estimate_richness(cag81_physeq, measures = "Shannon")!=0,1,0),
  "BMI" = phyloseq::sample_data(cag81_physeq)$BMI,
  "bmi_class" = phyloseq::sample_data(cag81_physeq)$bmi_class,
  "country" = phyloseq::sample_data(cag81_physeq)$country,
  "dataset_name" = phyloseq::sample_data(cag81_physeq)$dataset_name,
  "sex" = phyloseq::sample_data(cag81_physeq)$gender,
  "age" = phyloseq::sample_data(cag81_physeq)$age,
  "non_westernized" = phyloseq::sample_data(cag81_physeq)$non_westernized)
colnames(cag81_global_adiv) <- c("Observed", "shannon", "shannon_bin", "bmi", "bmi_class", "country", "dataset_name", "sex", "age", "non_westernized")

# BINARY MODEL (glm.bin)
Anova(glm(shannon_bin ~ bmi + sex + age, data = cag81_global_adiv, family = binomial(link = logit)))

# CONTINUOUS MODEL (glm.cont)
# remove columns where all values are 0 (for numeric values)   
cag81_shannon_cont <- cag81_global_adiv[cag81_global_adiv$shannon >0,]
Anova(glm(shannon ~ bmi + sex + age, data = cag81_shannon_cont, family = Gamma(link = log)))

# scatterplot
ggplot(cag81_shannon_cont, aes(x = bmi, y = shannon, color = bmi_class, shape = non_westernized)) +
  geom_point()


### g__COE1 (Lachnospiraceae)
# Calculate alpha diversity for COE1
coe1_global_adiv <- data.frame(
  "Observed" = phyloseq::estimate_richness(coe1_physeq, measures = "Observed"),
  "shannon" = phyloseq::estimate_richness(coe1_physeq, measures = "Shannon"),
  "shannon_bin" = ifelse(phyloseq::estimate_richness(coe1_physeq, measures = "Shannon")!=0,1,0),
  "BMI" = phyloseq::sample_data(coe1_physeq)$BMI,
  "bmi_class" = phyloseq::sample_data(coe1_physeq)$bmi_class,
  "country" = phyloseq::sample_data(coe1_physeq)$country,
  "dataset_name" = phyloseq::sample_data(coe1_physeq)$dataset_name,
  "sex" = phyloseq::sample_data(coe1_physeq)$gender,
  "age" = phyloseq::sample_data(coe1_physeq)$age,
  "non_westernized" = phyloseq::sample_data(coe1_physeq)$non_westernized)
colnames(coe1_global_adiv) <- c("Observed", "shannon", "shannon_bin", "bmi", "bmi_class", "country", "dataset_name", "sex", "age", "non_westernized")

# BINARY MODEL (glm.bin)
Anova(glm(shannon_bin ~ bmi + sex + age, data = coe1_global_adiv, family = binomial(link = logit)))

# CONTINUOUS MODEL (glm.cont)
# remove columns where all values are 0 (for numeric values)   
coe1_shannon_cont <- coe1_global_adiv[coe1_global_adiv$shannon >0,]
Anova(glm(shannon ~ bmi + sex + age, data = coe1_shannon_cont, family = Gamma(link = log)))

# scatterplot
ggplot(coe1_shannon_cont, aes(x = bmi, y = shannon, color = bmi_class, shape = non_westernized)) +
  geom_point()


### g__Ruminococcus_A (Lachnospiraceae)
# Calculate alpha diversity for Ruminococcus_A
ruminoA_global_adiv <- data.frame(
  "Observed" = phyloseq::estimate_richness(ruminoA_physeq, measures = "Observed"),
  "shannon" = phyloseq::estimate_richness(ruminoA_physeq, measures = "Shannon"),
  "shannon_bin" = ifelse(phyloseq::estimate_richness(ruminoA_physeq, measures = "Shannon")!=0,1,0),
  "BMI" = phyloseq::sample_data(ruminoA_physeq)$BMI,
  "bmi_class" = phyloseq::sample_data(ruminoA_physeq)$bmi_class,
  "country" = phyloseq::sample_data(ruminoA_physeq)$country,
  "dataset_name" = phyloseq::sample_data(ruminoA_physeq)$dataset_name,
  "sex" = phyloseq::sample_data(ruminoA_physeq)$gender,
  "age" = phyloseq::sample_data(ruminoA_physeq)$age,
  "non_westernized" = phyloseq::sample_data(ruminoA_physeq)$non_westernized)
colnames(ruminoA_global_adiv) <- c("Observed", "shannon", "shannon_bin", "bmi", "bmi_class", "country", "dataset_name", "sex", "age", "non_westernized")

# BINARY MODEL (glm.bin)
Anova(glm(shannon_bin ~ bmi + sex + age, data = ruminoA_global_adiv, family = binomial(link = logit)))

# CONTINUOUS MODEL (glm.cont)
# remove columns where all values are 0 (for numeric values)   
ruminoA_shannon_cont <- ruminoA_global_adiv[ruminoA_global_adiv$shannon >0,]
Anova(glm(shannon ~ bmi + sex + age, data = ruminoA_shannon_cont, family = Gamma(link = log)))

# scatterplot
ggplot(ruminoA_shannon_cont, aes(x = bmi, y = shannon, color = bmi_class, shape = non_westernized)) +
  geom_point()
