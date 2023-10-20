#source("../scripts/phylotags_preprocessing.R")

idlevel=100 #this matters for testing the results at different levels of identity (phylotags_suppl.R)
ssu_rare <- ssu_rare_100
dnaK_rare <- dnaK_rare_100
gyrB_rare <- gyrB_rare_100

######################
# ALPHA DIVERSITY ----
######################

# Generate a data.frame with alpha diversity measures
ssu_adiv <- data.frame(
  "Observed" = phyloseq::estimate_richness(ssu_rare, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ssu_rare, measures = "Shannon"),
  #  "PD" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(dnaK_rare)))), tree = phyloseq::phy_tree(dnaK_rare))[, 1],
  #  "CHS" = phyloseq::sample_data(dnaK_rare)$chs_class,
  "BMI" = phyloseq::sample_data(ssu_rare)$bmi,
  "bmi_class" = factor(sample_data(ssu_rare)$bmi_class,levels = c("Lean","Overweight","Obese")),
  "city" = phyloseq::sample_data(ssu_rare)$city,
  "sex" = phyloseq::sample_data(ssu_rare)$sex,
  "age_range" = phyloseq::sample_data(ssu_rare)$age_range,
  "stool_consistency" = phyloseq::sample_data(ssu_rare)$stool_consistency,
  "medicament" = phyloseq::sample_data(ssu_rare)$medicament)
# head(ssu_adiv)

dnaK_adiv <- data.frame(
  "Observed" = phyloseq::estimate_richness(dnaK_rare, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(dnaK_rare, measures = "Shannon"),
  #  "PD" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(dnaK_rare)))), tree = phyloseq::phy_tree(dnaK_rare))[, 1],
  #  "CHS" = phyloseq::sample_data(dnaK_rare)$chs_class,
  "BMI" = phyloseq::sample_data(dnaK_rare)$bmi,
  "bmi_class" = factor(sample_data(dnaK_rare)$bmi_class,levels = c("Lean","Overweight","Obese")),
  "city" = phyloseq::sample_data(dnaK_rare)$city,
  "sex" = phyloseq::sample_data(dnaK_rare)$sex,
  "age_range" = phyloseq::sample_data(dnaK_rare)$age_range,
  "stool_consistency" = phyloseq::sample_data(dnaK_rare)$stool_consistency,
  "medicament" = phyloseq::sample_data(dnaK_rare)$medicament)
# head(dnaK_adiv)

gyrB_adiv <- data.frame(
  "Observed" = phyloseq::estimate_richness(gyrB_rare, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(gyrB_rare, measures = "Shannon"),
  #  "PD" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(dnaK_rare)))), tree = phyloseq::phy_tree(dnaK_rare))[, 1],
  #  "CHS" = phyloseq::sample_data(dnaK_rare)$chs_class,
  "BMI" = phyloseq::sample_data(gyrB_rare)$bmi,
  "bmi_class" = factor(sample_data(gyrB_rare)$bmi_class,levels = c("Lean","Overweight","Obese")),
  "city" = phyloseq::sample_data(gyrB_rare)$city,
  "sex" = phyloseq::sample_data(gyrB_rare)$sex,
  "age_range" = phyloseq::sample_data(gyrB_rare)$age_range,
  "stool_consistency" = phyloseq::sample_data(gyrB_rare)$stool_consistency,
  "medicament" = phyloseq::sample_data(gyrB_rare)$medicament)
# head(gyrB_adiv)


### Some graphs and tests per covariates
# City of origin
pirateplot(formula = Shannon ~ city,
           data = ssu_adiv,
           pal = "google",
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           xlab = "City of origin",
           main = "16S rRNA")
kruskal.test(Shannon ~ city, data = ssu_adiv)

pirateplot(formula = Shannon ~ city,
           data = dnaK_adiv,
           pal = "google",
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           xlab = "City of origin",
           main = "dnaK")
kruskal.test(Shannon ~ city, data = dnaK_adiv)

pirateplot(formula = Shannon ~ city,
           data = gyrB_adiv,
           pal = "google",
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           xlab = "City of origin",
           main = "gyrB")
kruskal.test(Shannon ~ city, data = gyrB_adiv)

# Sex
pirateplot(formula = Shannon ~ sex,
           data = ssu_adiv,
           pal = "google",
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           xlab = "Sex at birth",
           main = "16S rRNA")
wilcox.test(Shannon ~ sex, data = ssu_adiv)

pirateplot(formula = Shannon ~ sex,
           data = dnaK_adiv,
           pal = "google",
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           xlab = "Sex at birth",
           main = "dnaK")
wilcox.test(Shannon ~ sex, data = dnaK_adiv)

pirateplot(formula = Shannon ~ sex,
           data = gyrB_adiv,
           pal = "google",
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           xlab = "Sex at birth",
           main = "gyrB")
wilcox.test(Shannon ~ sex, data = gyrB_adiv)

# Age range
pirateplot(formula = Shannon ~ age_range,
           data = ssu_adiv,
           pal = "google",
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           xlab = "Age category (years)",
           main = "16S rRNA")
wilcox.test(Shannon ~ age_range, data = ssu_adiv)

pirateplot(formula = Shannon ~ age_range,
           data = dnaK_adiv,
           pal = "google",
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           xlab = "Age category (years)",
           main = "dnaK")
wilcox.test(Shannon ~ age_range, data = dnaK_adiv)

pirateplot(formula = Shannon ~ age_range,
           data = gyrB_adiv,
           pal = "google",
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           xlab = "Age category (years)",
           main = "gyrB")
wilcox.test(Shannon ~ age_range, data = gyrB_adiv)

# Stool consistency
pirateplot(formula = Shannon ~ stool_consistency,
           data = ssu_adiv,
           pal = "google",
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           xlab = "Stool consistency",
           main = "16S rRNA")
kruskal.test(Shannon ~ stool_consistency, data = ssu_adiv)

pirateplot(formula = Shannon ~ stool_consistency,
           data = dnaK_adiv,
           pal = "google",
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           xlab = "Stool consistency",
           main = "dnaK")
kruskal.test(Shannon ~ stool_consistency, data = dnaK_adiv)

pirateplot(formula = Shannon ~ stool_consistency,
           data = gyrB_adiv,
           pal = "google",
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           xlab = "Stool consistency",
           main = "gyrB")
kruskal.test(Shannon ~ stool_consistency, data = gyrB_adiv)

# Medicament
pirateplot(formula = Shannon ~ medicament,
           data = ssu_adiv,
           pal = "google",
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           xlab = "Medication use",
           main = "16S rRNA")
wilcox.test(Shannon ~ medicament, data = ssu_adiv)

pirateplot(formula = Shannon ~ medicament,
           data = dnaK_adiv,
           pal = "google",
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           xlab = "Medication use",
           main = "dnaK")
wilcox.test(Shannon ~ medicament, data = dnaK_adiv)

pirateplot(formula = Shannon ~ medicament,
           data = gyrB_adiv,
           pal = "google",
           gl = 0, point.o = 0.8,
           inf.method = "ci",
           xlab = "Medication use",
           main = "gyrB")
wilcox.test(Shannon ~ medicament, data = gyrB_adiv)


###################################################
## SUMMARY OF SHANNON STATS FOR DIFFERENT TAXA ----
###################################################

### Function to get automated phyloseq with tax level and choice within
# https://github.com/joey711/phyloseq/issues/1471
fn <- function(ps, level, choice) {
  x <- paste0(level,"==","'",choice,"'")
  oldTax <- data.frame(tax_table(ps))
  newTax <- subset(oldTax, eval(parse(text=x)))
  tax_table(ps) <- tax_table(as.matrix(newTax))
  return(ps)
}

# Example
fn(dnaK_rare, level="family",choice="f__Acidaminococcaceae")

## FUNCTION MARKER_TAX_STATS
#this function takes marker and level of taxonomy as argument and returns
# 1) dataframe of Shannon values for each family
# 2) list of anovas 
# 3) list of glm 

#include idlevel clustering: 100, 97, 85 (for other than 100 run phylotags_suppl.R)

# Comment (#) and uncomment lines depending on the covariates you want to consider in your glm
marker_tax_stats = function(marker=NULL,clustering=NULL,level_tax=NULL){
  #take the rarefied phyloseq object for the marker
  phyobj=get(paste(marker,paste("rare",clustering,sep="_"),sep="_")) 
  #get a list of the families present
  taxa_list = levels(as.factor(tax_table(phyobj)[,level_tax]))
  #initialize a matrix with equal samples as the phyloseq object
  shannon_family.df = data.frame(matrix(nrow=ncol(otu_table(phyobj)),ncol=0))
  shannonbinary_family.df = data.frame(matrix(nrow=ncol(otu_table(phyobj)),ncol=0))
  observed_family.df = data.frame(matrix(nrow=ncol(otu_table(phyobj)),ncol=0))
  
  #loop over each taxon and calculate the Shannon index for each sample, then bind to a dataframe
  family_names = NULL
  for (tx in taxa_list){
    #  print(tx)
    shannon_family.df <- cbind(shannon_family.df,phyloseq::estimate_richness(fn(phyobj, level=level_tax,choice=tx),measures = "Shannon"))
    shannonbinary_family.df <- cbind(shannonbinary_family.df,ifelse(phyloseq::estimate_richness(fn(phyobj, level=level_tax,choice=tx),measures = "Shannon")!=0,1,0))
    observed_family.df <- cbind(observed_family.df,phyloseq::estimate_richness(fn(phyobj, level=level_tax,choice=tx),measures = "Observed"))
    
    # fix names that have stripe or space
    tx[grepl("-",tx)] <- str_replace_all(tx, "-","_")
    tx[grepl(" ", tx)] <- str_replace_all(tx, " ","_")
    tx[grepl("^[[:digit:]]+", tx)] <- paste("f",tx[grepl("^[[:digit:]]+", tx)],sep="")
    family_names <- append(family_names,tx)
  }
  
  # Fix names with non allowed characters    
  # fix names that have stripe or space
  # family_names[grepl("-",family_names)] <- str_replace_all(tx, "-","_")
  # spaces removed form names
  # family_names[grepl(" ", family_names)] <- family_names[grepl("", family_names)]
  #R does no deal well with strings starting with number, so here I add a letter to those
  # family_names[grepl("^[[:digit:]]+", family_names)] <- paste("f",family_names[grepl("^[[:digit:]]+", family_names)],sep="")
  
  colnames(shannon_family.df) <- family_names 
  colnames(shannonbinary_family.df) <- family_names
  colnames(observed_family.df) <- family_names
  
  #add BMI as explanatory variable
  shannon_family.df$bmi <- sample_data(phyobj)$bmi
  shannonbinary_family.df$bmi <- sample_data(phyobj)$bmi
  observed_family.df$bmi <- sample_data(phyobj)$bmi
  
#  shannon_family.df$bmi_class <- factor(sample_data(phyobj)$bmi_class,levels = c("Lean","Overweight","Obese"))
#  shannonbinary_family.df$bmi_class <- factor(sample_data(phyobj)$bmi_class,levels = c("Lean","Overweight","Obese"))
#  observed_family.df$bmi_class <- factor(sample_data(phyobj)$bmi_class,levels = c("Lean","Overweight","Obese"))
  
  #add city as explanatory variable
  shannon_family.df$city <- factor(sample_data(phyobj)$city)
  shannonbinary_family.df$city <- factor(sample_data(phyobj)$city)
  observed_family.df$city <- factor(sample_data(phyobj)$city)
  
  #add sex as explanatory variable
#  shannon_family.df$sex <- factor(sample_data(phyobj)$sex)
#  shannonbinary_family.df$sex <- factor(sample_data(phyobj)$sex)
#  observed_family.df$sex <- factor(sample_data(phyobj)$sex)
  
  #add age_range as explanatory variable
#  shannon_family.df$age_range <- factor(sample_data(phyobj)$age_range)
#  shannonbinary_family.df$age_range <- factor(sample_data(phyobj)$age_range)
#  observed_family.df$age_range <- factor(sample_data(phyobj)$age_range)
  
  #add stool_consistency as explanatory variable
#  shannon_family.df$stool_consistency <- factor(sample_data(phyobj)$stool_consistency)
#  shannonbinary_family.df$stool_consistency <- factor(sample_data(phyobj)$stool_consistency)
#  observed_family.df$stool_consistency <- factor(sample_data(phyobj)$stool_consistency)
  
  #add medicament as explanatory variable
#  shannon_family.df$medicament <- factor(sample_data(phyobj)$medicament)
#  shannonbinary_family.df$medicament <- factor(sample_data(phyobj)$medicament)
#  observed_family.df$medicament <- factor(sample_data(phyobj)$medicament)
  
  diversity.df <- list(shannon=shannon_family.df,shannon_bin=shannonbinary_family.df,observed=observed_family.df)
  return(diversity.df)
}


###############################################################
# Gamma-hurdle (zero and non-zero modeled "independently") ----
###############################################################

#https://seananderson.ca/2014/05/18/gamma-hurdle/
#two models that give
# 1) m1: the probability of non-zero 
# 2) given that is non-zero the response of y
# 1) binary model
# adjusted to city variable

### BINARY (glm.bin) ---
## 16S rRNA
ssu_shannon_bin.fam <- marker_tax_stats("ssu",idlevel,"family")$shannon_bin
ssu_shannon_bin.gen <- marker_tax_stats("ssu",idlevel,"genus")$shannon_bin

# include more or less covariates in this list depending on your models (e.g., bmi, city, sex, age_range, etc.)
ssu_vars.fam = names(ssu_shannon_bin.fam)[!names(ssu_shannon_bin.fam) %in% c("bmi","city")]
ssu_vars.gen = names(ssu_shannon_bin.gen)[!names(ssu_shannon_bin.gen) %in% c("bmi","city")]

ssu_glmbin_fam.df <- lapply(setNames(ssu_vars.fam, ssu_vars.fam),function(s) {
  tryCatch(glm(as.formula(paste(as.character(s), " ~ bmi + city")),
               data = ssu_shannon_bin.fam, family = binomial(link = logit)),
           error=function(e) NULL)})

ssu_glmbin_gen.df <- lapply(setNames(ssu_vars.gen, ssu_vars.gen),function(s) {
  tryCatch(glm(as.formula(paste(as.character(s), " ~ bmi + city")),
               data = ssu_shannon_bin.gen, family = binomial(link = logit)), 
           error=function(e) NULL)})
  
ssu_glmbin_fam.mapdf <- map_df(ssu_glmbin_fam.df, tidy, .id="taxa")
ssu_glmbin_gen.mapdf <- map_df(ssu_glmbin_gen.df, tidy, .id="taxa")

ssu_glmbin_fam.sigtab <- ssu_glmbin_fam.mapdf %>%
  filter(term=="bmi")
ssu_glmbin_fam.sigtab$p.adjust <-p.adjust(ssu_glmbin_fam.sigtab$p.value, "fdr")
ssu_glmbin_fam.sigtab %>%
  filter(p.value<0.05)

ssu_glmbin_gen.sigtab <- ssu_glmbin_gen.mapdf %>%
  filter(term=="bmi")
ssu_glmbin_gen.sigtab$p.adjust <-p.adjust(ssu_glmbin_gen.sigtab$p.value, "fdr")
ssu_glmbin_gen.sigtab %>%
  filter(p.value<0.05)


## dnaK
dnaK_shannon_bin.fam <- marker_tax_stats("dnaK",idlevel,"family")$shannon_bin
dnaK_shannon_bin.gen <- marker_tax_stats("dnaK",idlevel,"genus")$shannon_bin

# include more or less covariates in this list depending on your models (e.g., bmi, city, sex, age_range, etc.)
dnaK_vars.fam = names(dnaK_shannon_bin.fam)[!names(dnaK_shannon_bin.fam) %in% c("bmi","city")]
dnaK_vars.gen = names(dnaK_shannon_bin.gen)[!names(dnaK_shannon_bin.gen) %in% c("bmi","city")]

dnaK_glmbin_fam.df <- lapply(setNames(dnaK_vars.fam, dnaK_vars.fam),function(s) {
  tryCatch(glm(as.formula(paste(as.character(s), " ~ bmi + city")),
               data = dnaK_shannon_bin.fam, family = binomial(link = logit)), 
           error=function(e) NULL)})

dnaK_glmbin_gen.df <- lapply(setNames(dnaK_vars.gen, dnaK_vars.gen),function(s) {
  tryCatch(glm(as.formula(paste(as.character(s), " ~ bmi + city")),
               data = dnaK_shannon_bin.gen, family = binomial(link = logit)), 
           error=function(e) NULL)})

dnaK_glmbin_fam.mapdf <- map_df(dnaK_glmbin_fam.df, tidy, .id="taxa")
dnaK_glmbin_gen.mapdf <- map_df(dnaK_glmbin_gen.df, tidy, .id="taxa")

# Adjust p-values within the families covered by the primers
families_dnaK <- levels(as.factor(dnaK_taxmat_dada2[which(dnaK_taxmat_dada2[,"family"] %in% "f__Oscillospiraceae" | dnaK_taxmat_dada2[,"family"] %in% "f__Ruminococcaceae" | dnaK_taxmat_dada2[,"family"] %in% "f__Acutalibacteraceae"),"family"]))
families_dnaK <- str_replace_all(families_dnaK, "-","_")
dnaK_glmbin_fam.families_dnaK <- dnaK_glmbin_fam.mapdf%>% 
  filter(taxa %in% families_dnaK & term=="bmi") %>%
  mutate(padj=p.adjust(p.value,method = "fdr"))
dnaK_glmbin_fam.families_dnaK %>%
  filter(p.value<0.05)

oscillo_dnaK <- levels(as.factor(dnaK_taxmat_dada2[which(dnaK_taxmat_dada2[,"family"] %in% "f__Oscillospiraceae"),"genus"]))
oscillo_dnaK <- str_replace_all(oscillo_dnaK, "-","_")
dnaK_glmbin_gen.oscillo_dnaK <- dnaK_glmbin_gen.mapdf%>% 
  filter(taxa %in% oscillo_dnaK & term=="bmi") %>%
  mutate(padj=p.adjust(p.value,method = "fdr"))
dnaK_glmbin_gen.oscillo_dnaK %>%
  filter(p.value<0.05)

rumino_dnaK <- levels(as.factor(dnaK_taxmat_dada2[which(dnaK_taxmat_dada2[,"family"] %in% "f__Ruminococcaceae"),"genus"]))
rumino_dnaK <- str_replace_all(rumino_dnaK, "-","_")
dnaK_glmbin_gen.rumino_dnaK <- dnaK_glmbin_gen.mapdf%>% 
  filter(taxa %in% rumino_dnaK & term=="bmi") %>%
  mutate(padj=p.adjust(p.value,method = "fdr"))
dnaK_glmbin_gen.rumino_dnaK %>%
  filter(p.value<0.05)

acutali_dnaK <- levels(as.factor(dnaK_taxmat_dada2[which(dnaK_taxmat_dada2[,"family"] %in% "f__Acutalibacteraceae"),"genus"]))
acutali_dnaK <- str_replace_all(acutali_dnaK, "-","_")
dnaK_glmbin_gen.acutali_dnaK <- dnaK_glmbin_gen.mapdf%>% 
  filter(taxa %in% acutali_dnaK & term=="bmi") %>%
  mutate(padj=p.adjust(p.value,method = "fdr"))
dnaK_glmbin_gen.acutali_dnaK %>%
  filter(p.value<0.05)


## gyrB
gyrB_shannon_bin.fam <- marker_tax_stats("gyrB",idlevel,"family")$shannon_bin
gyrB_shannon_bin.gen <- marker_tax_stats("gyrB",idlevel,"genus")$shannon_bin
      
# include more or less covariates in this list depending on your models (e.g., bmi, city, sex, age_range, etc.)
gyrB_vars.fam = names(gyrB_shannon_bin.fam)[!names(gyrB_shannon_bin.fam) %in% c("bmi","city")]
gyrB_vars.gen = names(gyrB_shannon_bin.gen)[!names(gyrB_shannon_bin.gen) %in% c("bmi","city")]

gyrB_glmbin_fam.df <- lapply(setNames(gyrB_vars.fam, gyrB_vars.fam),function(s) {
  tryCatch(glm(as.formula(paste(as.character(s), " ~ bmi + city")),
               data = gyrB_shannon_bin.fam, family = binomial(link = logit)), 
           error=function(e) NULL)})

gyrB_glmbin_gen.df <- lapply(setNames(gyrB_vars.gen, gyrB_vars.gen),function(s) {
  tryCatch(glm(as.formula(paste(as.character(s), " ~ bmi + city")),
               data = gyrB_shannon_bin.gen, family = binomial(link = logit)), 
           error=function(e) NULL)})

gyrB_glmbin_fam.mapdf <- map_df(gyrB_glmbin_fam.df, tidy, .id="taxa")
gyrB_glmbin_gen.mapdf <- map_df(gyrB_glmbin_gen.df, tidy, .id="taxa")
  

# Adjust p-values within the families covered by the primers
families_gyrB <- levels(as.factor(gyrB_taxmat_dada2[which(gyrB_taxmat_dada2[,"family"] %in% "f__Lachnospiraceae"),"family"]))
families_gyrB <- str_replace_all(families_gyrB, "-","_")
gyrB_glmbin_fam.families_gyrB <- gyrB_glmbin_fam.mapdf%>% 
  filter(taxa %in% families_gyrB & term=="bmi") %>%
  mutate(padj=p.adjust(p.value,method = "fdr"))
gyrB_glmbin_fam.families_gyrB %>%
  filter(p.value<0.05)

lachno_gyrB <- levels(as.factor(gyrB_taxmat_dada2[which(gyrB_taxmat_dada2[,"family"] %in% "f__Lachnospiraceae"),"genus"]))
lachno_gyrB <- str_replace_all(lachno_gyrB, "-","_")
gyrB_glmbin_gen.lachno_gyrB <- gyrB_glmbin_gen.mapdf%>% 
  filter(taxa %in% lachno_gyrB & term=="bmi") %>%
  mutate(padj=p.adjust(p.value,method = "fdr"))
gyrB_glmbin_gen.lachno_gyrB %>%
  filter(p.value<0.05)



## CONTINOUS (glm.cont) ----
## 16S rRNA
ssu_shannon.fam <- marker_tax_stats("ssu",idlevel,"family")$shannon
ssu_shannon.gen <- marker_tax_stats("ssu",idlevel,"genus")$shannon

# first remove columns where all values are 0 (for numeric values)   
ssu_shannon_cont.fam <- ssu_shannon.fam %>%
  select_if(~ !is.numeric(.) || sum(.) != 0)
ssu_shannon_cont.gen <- ssu_shannon.gen %>%
  select_if(~ !is.numeric(.) || sum(.) != 0)

# include more or less covariates in this list depending on your models (e.g., bmi, city, sex, age_range, etc.)
ssu_vars.fam = names(ssu_shannon_cont.fam)[!names(ssu_shannon_cont.fam) %in% c("bmi","city")]
ssu_vars.gen = names(ssu_shannon_cont.gen)[!names(ssu_shannon_cont.gen) %in% c("bmi","city")]

# second drop also variables with less than 3 values (at least two values different to 0), otherwise there will be errors when attempting
# regression model using a predictor variable with 2 or less values
ssu_keep.fam <- names(which(sapply(lapply(ssu_shannon_cont.fam[,ssu_vars.fam], unique), length) > 3))
ssu_keep.gen <- names(which(sapply(lapply(ssu_shannon_cont.gen[,ssu_vars.gen], unique), length) > 3))

ssu_glmcont_fam.df <- lapply(setNames(ssu_keep.fam, ssu_keep.fam),function(s) {
  tryCatch(glm(as.formula(paste(as.character(s), " ~ bmi + city")),
               data = ssu_shannon_cont.fam[ssu_shannon_cont.fam[,s]!=0,], family = Gamma(link = log)),
           error=function(e) NULL)})

ssu_glmcont_gen.df <- lapply(setNames(ssu_keep.gen, ssu_keep.gen),function(s) {
  tryCatch(glm(as.formula(paste(as.character(s), " ~ bmi + city")),
               data = ssu_shannon_cont.gen[ssu_shannon_cont.gen[which(colnames(ssu_shannon_cont.gen)==s)]!=0,], family = Gamma(link = log)),
           error=function(e) NULL)})

ssu_glmcont_fam.mapdf <- map_df(ssu_glmcont_fam.df, tidy, .id="taxa")
ssu_glmcont_gen.mapdf <- map_df(ssu_glmcont_gen.df, tidy, .id="taxa")

ssu_glmcont_fam.sigtab <- ssu_glmcont_fam.mapdf %>%
  filter(term=="bmi")
ssu_glmcont_fam.sigtab$p.adjust <-p.adjust(ssu_glmcont_fam.sigtab$p.value, "fdr")
ssu_glmcont_fam.sigtab %>%
  filter(p.value<0.055)

ssu_glmcont_gen.sigtab <- ssu_glmcont_gen.mapdf %>%
  filter(term=="bmi")
ssu_glmcont_gen.sigtab$p.adjust <-p.adjust(ssu_glmcont_gen.sigtab$p.value, "fdr")
ssu_glmcont_gen.sigtab %>%
  filter(p.value<0.055)


## dnaK 
dnaK_shannon.fam <- marker_tax_stats("dnaK",idlevel,"family")$shannon
dnaK_shannon.gen <- marker_tax_stats("dnaK",idlevel,"genus")$shannon

# first remove columns where all values are 0 (for numeric values)
dnaK_shannon_cont.fam <- dnaK_shannon.fam %>%
  select_if(~ !is.numeric(.) || sum(.) != 0)
dnaK_shannon_cont.gen <- dnaK_shannon.gen %>%
  select_if(~ !is.numeric(.) || sum(.) != 0)

# include more or less covariates in this list depending on your models (e.g., bmi, city, sex, age_range, etc.)
dnaK_vars.fam = names(dnaK_shannon_cont.fam)[!names(dnaK_shannon_cont.fam) %in% c("bmi","city")]
dnaK_vars.gen = names(dnaK_shannon_cont.gen)[!names(dnaK_shannon_cont.gen) %in% c("bmi","city")]

# second drop also variables with less than 3 values (at least two values different to 0), otherwise there will be errors when attempti$
# regression model using a predictor variable with 2 or less values
dnaK_keep.fam <- names(which(sapply(lapply(dnaK_shannon_cont.fam[,dnaK_vars.fam], unique), length) > 3))
dnaK_keep.gen <- names(which(sapply(lapply(dnaK_shannon_cont.gen[,dnaK_vars.gen], unique), length) > 3))

dnaK_glmcont_fam.df <- lapply(setNames(dnaK_keep.fam, dnaK_keep.fam),function(s) {
  tryCatch(glm(as.formula(paste(as.character(s), " ~ bmi + city")),
               data = dnaK_shannon_cont.fam[dnaK_shannon_cont.fam[,s]!=0,], family = Gamma(link = log)),
           error=function(e) NULL)})

dnaK_glmcont_gen.df <- lapply(setNames(dnaK_keep.gen, dnaK_keep.gen),function(s) {
  tryCatch(glm(as.formula(paste(as.character(s), " ~ bmi + city")),
               data = dnaK_shannon_cont.gen[dnaK_shannon_cont.gen[,s]!=0,], family = Gamma(link = log)), 
           error=function(e) NULL)})
  
dnaK_glmcont_fam.mapdf <- map_df(dnaK_glmcont_fam.df, tidy, .id="taxa")
dnaK_glmcont_gen.mapdf <- map_df(dnaK_glmcont_gen.df, tidy, .id="taxa")

# Adjust p-values within the families covered by the primers
families_dnaK <- levels(as.factor(dnaK_taxmat_dada2[which(dnaK_taxmat_dada2[,"family"] %in% "f__Oscillospiraceae" | dnaK_taxmat_dada2[,"family"] %in% "f__Ruminococcaceae" | dnaK_taxmat_dada2[,"family"] %in% "f__Acutalibacteraceae"),"family"]))
families_dnaK <- str_replace_all(families_dnaK, "-","_")
dnaK_glmcont_fam.families_dnaK <- dnaK_glmcont_fam.mapdf%>% 
  filter(taxa %in% families_dnaK & term=="bmi") %>%
  mutate(padj=p.adjust(p.value,method = "fdr"))
dnaK_glmcont_fam.families_dnaK %>%
  filter(p.value<0.055)

oscillo_dnaK <- levels(as.factor(dnaK_taxmat_dada2[which(dnaK_taxmat_dada2[,"family"] %in% "f__Oscillospiraceae"),"genus"]))
oscillo_dnaK <- str_replace_all(oscillo_dnaK, "-","_")
dnaK_glmcont_gen.oscillo_dnaK <- dnaK_glmcont_gen.mapdf%>% 
  filter(taxa %in% oscillo_dnaK & term=="bmi") %>%
  mutate(padj=p.adjust(p.value,method = "fdr"))
dnaK_glmcont_gen.oscillo_dnaK %>%
  filter(p.value<0.055)

rumino_dnaK <- levels(as.factor(dnaK_taxmat_dada2[which(dnaK_taxmat_dada2[,"family"] %in% "f__Ruminococcaceae"),"genus"]))
rumino_dnaK <- str_replace_all(rumino_dnaK, "-","_")
dnaK_glmcont_gen.rumino_dnaK <- dnaK_glmcont_gen.mapdf%>% 
  filter(taxa %in% rumino_dnaK & term=="bmi") %>%
  mutate(padj=p.adjust(p.value,method = "fdr"))
dnaK_glmcont_gen.rumino_dnaK %>%
  filter(p.value<0.055)

acutali_dnaK <- levels(as.factor(dnaK_taxmat_dada2[which(dnaK_taxmat_dada2[,"family"] %in% "f__Acutalibacteraceae"),"genus"]))
acutali_dnaK <- str_replace_all(acutali_dnaK, "-","_")
dnaK_glmcont_gen.acutali_dnaK <- dnaK_glmcont_gen.mapdf%>% 
  filter(taxa %in% acutali_dnaK & term=="bmi") %>%
  mutate(padj=p.adjust(p.value,method = "fdr"))
dnaK_glmcont_gen.acutali_dnaK %>%
  filter(p.value<0.055)


## gyrB
gyrB_shannon.fam <- marker_tax_stats("gyrB",idlevel,"family")$shannon
gyrB_shannon.gen <- marker_tax_stats("gyrB",idlevel,"genus")$shannon

# first remove columns where all values are 0 (for numeric values)
gyrB_shannon_cont.fam <- gyrB_shannon.fam %>%
  select_if(~ !is.numeric(.) || sum(.) != 0)
gyrB_shannon_cont.gen <- gyrB_shannon.gen %>%
  select_if(~ !is.numeric(.) || sum(.) != 0)
      
# include more or less covariates in this list depending on your models (e.g., bmi, city, sex, age_range, etc.)
gyrB_vars.fam = names(gyrB_shannon_cont.fam)[!names(gyrB_shannon_cont.fam) %in% c("bmi","city")]
gyrB_vars.gen = names(gyrB_shannon_cont.gen)[!names(gyrB_shannon_cont.gen) %in% c("bmi","city")]

# second drop also variables with less than 3 values (at least two values different to 0), otherwise there will be errors when attempti$
# regression model using a predictor variable with 2 or less values
gyrB_keep.fam <- names(which(sapply(lapply(gyrB_shannon_cont.fam[,gyrB_vars.fam], unique), length) > 3))
gyrB_keep.gen <- names(which(sapply(lapply(gyrB_shannon_cont.gen[,gyrB_vars.gen], unique), length) > 3))
  
gyrB_glmcont_fam.df <- lapply(setNames(gyrB_keep.fam, gyrB_keep.fam),function(s) {
  tryCatch(glm(as.formula(paste(as.character(s), " ~ bmi + city")),
               data = gyrB_shannon_cont.fam[gyrB_shannon_cont.fam[,s]!=0,], family = Gamma(link = log)), 
           error=function(e) NULL)})

gyrB_glmcont_gen.df <- lapply(setNames(gyrB_keep.gen, gyrB_keep.gen),function(s) {
  tryCatch(glm(as.formula(paste(as.character(s), " ~ bmi + city")),
               data = gyrB_shannon_cont.gen[gyrB_shannon_cont.gen[,s]!=0,], family = Gamma(link = log)),
           error=function(e) NULL)})

gyrB_glmcont_fam.mapdf <- map_df(gyrB_glmcont_fam.df, tidy, .id="taxa")
gyrB_glmcont_gen.mapdf <- map_df(gyrB_glmcont_gen.df, tidy, .id="taxa")

# Adjust p-values within the families covered by the primers
families_gyrB <- levels(as.factor(gyrB_taxmat_dada2[which(gyrB_taxmat_dada2[,"family"] %in% "f__Lachnospiraceae"),"family"]))
families_gyrB <- str_replace_all(families_gyrB, "-","_")
gyrB_glmcont_fam.families_gyrB <- gyrB_glmcont_fam.mapdf%>% 
  filter(taxa %in% families_gyrB & term=="bmi") %>%
  mutate(padj=p.adjust(p.value,method = "fdr"))
gyrB_glmcont_fam.families_gyrB %>%
  filter(p.value<0.055)

lachno_gyrB <- levels(as.factor(gyrB_taxmat_dada2[which(gyrB_taxmat_dada2[,"family"] %in% "f__Lachnospiraceae"),"genus"]))
lachno_gyrB <- str_replace_all(lachno_gyrB, "-","_")
gyrB_glmcont_gen.lachno_gyrB <- gyrB_glmcont_gen.mapdf%>% 
  filter(taxa %in% lachno_gyrB & term=="bmi") %>%
  mutate(padj=p.adjust(p.value,method = "fdr"))
gyrB_glmcont_gen.lachno_gyrB %>%
  filter(p.value<0.055)


#####################################
# Figure 4: GLMs diversity ~ bmi ----
#####################################

# Color blind friendly palette
okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

### Oscillospiraceae
# dnaK
oscillo_dnaK_fam.df <- cbind(dnaK_shannon_cont.fam %>% select(bmi,shannon_cont = "f__Oscillospiraceae"),
                             dnaK_shannon_bin.fam %>% select(shannon_bin = "f__Oscillospiraceae"))
oscillo_dnaK_fam.df$labelling = ifelse(oscillo_dnaK_fam.df$shannon_bin==1, ">0",oscillo_dnaK_fam.df$shannon_bin)

oscillo_dnaK_pval <- as.numeric(dnaK_glmcont_fam.families_dnaK[dnaK_glmcont_fam.families_dnaK["taxa"]=="f__Oscillospiraceae","p.value"])

# 16S rRNA
oscillo_ssu_fam.df <- cbind(ssu_shannon_cont.fam %>% select(bmi,shannon_cont = "Oscillospiraceae"),
                            ssu_shannon_bin.fam %>% select(shannon_bin = "Oscillospiraceae"))
oscillo_ssu_fam.df$labelling = ifelse(oscillo_ssu_fam.df$shannon_bin==1, ">0",oscillo_ssu_fam.df$shannon_bin)

oscillo_ssu_pval <- as.numeric(ssu_glmcont_fam.sigtab[ssu_glmcont_fam.sigtab["taxa"]=="Oscillospiraceae","p.value"])

# Oscillospiraceae dnaK + 16S rRNA in one plot
oscillo_plot <- plot_grid(ggplot(oscillo_ssu_fam.df, aes(bmi,shannon_cont,color=as.factor(labelling) )) + 
                            geom_point(size=2.5) +
                            theme(legend.title=element_blank()) +
                            theme(legend.position = "none")+
                            scale_colour_manual(name="shannon",values=okabe) +
                            ylim(c(0,4)) + 
                            xlim(c(18,40)) + 
                            geom_smooth( data=subset(oscillo_ssu_fam.df,shannon_bin==1),
                                         aes(bmi,shannon_cont),
                                         method = "glm", method.args = list(family = "Gamma"),
                                         se = FALSE,  linewidth=.75) + 
                            #annotate("text", x = 25, y = 4, label = paste("cont.pval",round(oscillo_ssu_pval,digits=4),sep="=")) +
                            ggtitle("Oscillospiraceae 16S rRNA") +
                            labs(x = bquote("BMI ("~kg/m^2 ~")"), y = "Shannon\n"),
                          
                          ggplot(oscillo_dnaK_fam.df, aes(bmi,shannon_cont,color=as.factor(labelling) )) + 
                            geom_point(size=2.5) +
                            theme(legend.title=element_blank()) +
                            theme(legend.position = "none")+
                            scale_colour_manual(name="shannon",values=okabe) +
                            ylim(c(0,4)) +
                            xlim(c(18,40)) +
                            geom_smooth(data=subset(oscillo_dnaK_fam.df,shannon_bin==1),
                                        aes(bmi,shannon_cont),
                                        method = "glm",method.args = list(family = "Gamma"),
                                        se = FALSE,  linewidth=.75) +  
                            #annotate("text", x = 25, y = 4,  label = paste("cont.pval",round(oscillo_dnaK_pval,digits=4),sep="=")) + 
                            ggtitle("Oscillospiraceae dnaK") + labs(x = bquote("BMI ("~kg/m^2 ~")"), y = "Shannon\n"),
                          nrow = 1, ncol = 2, label_size = 10)


### Lachnospiraceae
# gyrB
lachno_gyrB_fam.df <- cbind(gyrB_shannon_cont.fam %>% select(bmi,shannon_cont = "f__Lachnospiraceae"),
                            gyrB_shannon_bin.fam %>% select(shannon_bin = "f__Lachnospiraceae"))
lachno_gyrB_fam.df$labelling = ifelse(lachno_gyrB_fam.df$shannon_bin==1, ">0",lachno_gyrB_fam.df$shannon_bin)
lachno_gyrB_pval <- as.numeric(gyrB_glmcont_fam.families_gyrB[gyrB_glmcont_fam.families_gyrB["taxa"]=="f__Lachnospiraceae","p.value"])

# 16S rRNA
lachno_ssu_fam.df <- cbind(ssu_shannon_cont.fam %>% select(bmi,shannon_cont = "Lachnospiraceae"),
                           ssu_shannon_bin.fam %>% select(shannon_bin = "Lachnospiraceae"))
lachno_ssu_fam.df$labelling = ifelse(lachno_ssu_fam.df$shannon_bin==1, ">0",lachno_ssu_fam.df$shannon_bin)
lachno_ssu_pval <- as.numeric(ssu_glmcont_fam.sigtab[ssu_glmcont_fam.sigtab["taxa"]=="Lachnospiraceae","p.value"])#0.00365 #gyrB_glmcont_city_fam.sigtab$p.value

lachno_title <- ggdraw() + draw_label("Lachnospiraceae") #, fontface='bold')
lachno_plot <- plot_grid(ggplot(lachno_ssu_fam.df, aes(bmi,shannon_cont,color=as.factor(labelling) )) + 
                           geom_point(size=2.5) +
                           theme(legend.title=element_blank()) +
                           theme(legend.position = "none") +
                           scale_colour_manual(name="shannon",values=okabe) +
                           ylim(c(0,4)) + 
                           xlim(c(18,40)) + 
                           geom_smooth(data=subset(lachno_ssu_fam.df,shannon_bin==1), 
                                       aes(bmi,shannon_cont), 
                                       method = "glm", method.args = list(family = "Gamma"),
                                       se = FALSE, linewidth =.75) + 
                           #annotate("text", x = 25, y = 3.9, label = paste("cont.pval",round(lachno_ssu_pval,digits = 4),sep="=")) +  #"0.023"
                           ggtitle("Lachnospiraceae 16S rRNA") +
                           labs(x = bquote("BMI ("~kg/m^2 ~")"), y = "Shannon\n"),
                         ggplot(lachno_gyrB_fam.df, aes(bmi,shannon_cont,color=as.factor(labelling) )) + 
                           geom_point(size=2.5) +
                           theme(legend.title=element_blank()) +
                           theme(legend.position = "none")+
                           scale_colour_manual(name="shannon",values=okabe) +
                           ylim(c(0,4)) +
                           xlim(c(18,40)) +
                           geom_smooth(data=subset(lachno_gyrB_fam.df,shannon_bin==1),
                                       aes(bmi,shannon_cont), 
                                       method = "glm", method.args = list(family = "Gamma"),
                                       se = FALSE, linewidth =.75) + 
                           #annotate("text", x = 25, y = 3.9,  label = paste("cont.pval",round(lachno_gyrB_pval,digits = 4),sep="=")) + 
                           ggtitle("Lachnospiraceae gyrB") + labs(x = bquote("BMI ("~kg/m^2 ~")"), y = "Shannon\n"), 
                         nrow = 1, ncol = 2, label_size = 10)


### g__CAG-83 (Oscillospiraceae)
cag83 <- dnaK_shannon.gen %>% 
  select(bmi,shannon_cont = "g__CAG_83") %>% 
  mutate(shannon_bin = ifelse(shannon_cont==0,0,1)) %>%
  mutate(labelling = ifelse(shannon_cont==0,"0",">0"))

cag83_bin.boxplot <-ggplot(cag83)+
  geom_boxplot(size=1,aes(shannon_bin,bmi, color= labelling)) +
  #  geom_point(size=2.5,aes(bmi, shannon_bin , color= labelling)) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(name="shannon",values=okabe) +
  scale_x_discrete(name="Shannon", limits=c(0,1))+
  theme(legend.title=element_blank()) +
  theme(legend.position = "none") +
  #annotate("text", x = 0, y = 20,  label = paste("bin.pval",round(cag83.pval,digits = 4),sep="=")) +
  #annotate("text", x = 1.5, y = 20,  label = paste("bin.fdr",round(cag83.fdr,digits = 4),sep="=")) +
  ylim(c(18,40)) +
  labs(x = "Shannon \n", y = bquote("BMI ("~kg/m^2 ~")")) +
  ggtitle("g__CAG-83 dnaK")


### g__CAG-170 (Oscillospiraceae)
cag170 <- dnaK_shannon.gen %>% 
  select(bmi,shannon_cont = "g__CAG_170") %>% 
  mutate(shannon_bin = ifelse(shannon_cont==0,0,1)) %>%
  mutate(labelling = ifelse(shannon_cont==0,"0",">0"))

cag170_bin.boxplot <-ggplot(cag170)+
  geom_boxplot(size=1,aes(shannon_bin,bmi, color= labelling)) +
  #  geom_point(size=2.5,aes(bmi, shannon_bin , color= labelling)) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(name="shannon",values=okabe) +
  scale_x_discrete(name="Shannon", limits=c(0,1))+
  theme(legend.title=element_blank()) +
  theme(legend.position = "none") +
  #annotate("text", x = 0, y = 20,  label = paste("bin.pval",round(cag170.pval,digits = 4),sep="=")) +
  #annotate("text", x = 1.5, y = 20,  label = paste("bin.fdr",round(cag170.fdr,digits = 4),sep="=")) +
  ylim(c(18,40)) +
  labs(x = "Shannon \n", y = bquote("BMI ("~kg/m^2 ~")"))+
  ggtitle("g__CAG-170 dnaK")


Fig4 <- plot_grid(oscillo_plot, lachno_plot, plot_grid(cag83_bin.boxplot, cag170_bin.boxplot, ncol=2, nrow=1), ncol=1, nrow=3, labels = "AUTO")
