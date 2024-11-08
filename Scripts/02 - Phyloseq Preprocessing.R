# Author: Sunni Patton
# Last edited: 11/08/2024
# Title: Phyloseq preprocessing
# Overview: Preparing initial phyloseq object(s)

## Set seed ====
set.seed(123)

## Load libraries ====
library(here)
library(phyloseq)
library(vegan)

## Load data ====
read.csv(here::here("Output/00 - Read Preprocessing Output/metadata.csv")) -> metadata
readRDS(here::here("Output/00 - Read Preprocessing Output/seq_table_nobact.rds")) -> seq_table_nobact
readRDS(here::here("Output/00 - Read Preprocessing Output/taxa_nobact.rds")) -> taxa_nobact

## Make sample numbers into sample names to work with phyloseq ====
x <- metadata$SampleID
rownames(metadata) <- x

## Create phyloseq object ====
ps.All <- phyloseq(otu_table(seq_table_nobact, taxa_are_rows=FALSE), sample_data(metadata), tax_table(taxa_nobact))
saveRDS(ps.All, here::here("Output/02 - Phyloseq Preprocessing Output/ps.All.rds"))

## Make rarefaction curve ====
as.matrix(as.data.frame(ps.All@otu_table)) -> data
sort(sample_sums(ps.All))

# Make rarecurve to see how much diversity is captured if we rarefy to lowest sample depth (1084)
vegan::rarecurve(data, step = 100, label = FALSE, abline(v = 1084), col = "blue") # Essentially no diversity is captured
# Make rarecurve to see how much diversity is captured if we rarefy to 14205
vegan::rarecurve(data, step = 100, label = FALSE, abline(v = 14205, col = "red"), col = "black") # Most diversity is captured

## If we choose to rarefy to 14205, we automatically lose 2 samples and we still won't capture all the diversity 

## Remove negative controls and samples with ~1000 reads ====
subset_samples(ps.All, SampleID != "NC_1") -> ps.All
subset_samples(ps.All, SampleID != "NC_2") -> ps.All
subset_samples(ps.All, SampleID != "NC_3") -> ps.All
subset_samples(ps.All, SampleID != "NC_4") -> ps.All

subset_samples(ps.All, SampleID != "1_BAK_ACR_57_M19") -> ps.All
subset_samples(ps.All, SampleID != "1_BAK_ACR_57_A18") -> ps.All

# Save phyloseq object
saveRDS(ps.All, here::here("Output/02 - Phyloseq Preprocessing Output/ps.All.rds"))

## Run rrarefy() ====
# Get otu_table dataframe 
otu_table <- as.matrix(as.data.frame(ps.All@otu_table))
View(otu_table)

# Rarefy
rarefied_df <- rrarefy(otu_table, sample = 14205)
sort(rowSums(rarefied_df)) # All samples now have 14205 reads

# Make new phyloseq object after subsampling to even depth 
rare_samData <- ps.All@sam_data
rare_taxTable <- ps.All@tax_table
rare_otuTable <- rarefied_df

ps.rare <- phyloseq(otu_table(rare_otuTable, taxa_are_rows=FALSE),sample_data(rare_samData),tax_table(rare_taxTable))

# Make sure there aren't any taxa present that aren't actually in any sample
ps.rare <- prune_taxa(taxa_sums(ps.rare@otu_table) > 0, ps.rare) # Removed 201 taxa

summary(taxa_sums(ps.rare)) # minimum times a taxa appears: 1, 1st quartile: 3, mediamn: 20, mean: 1157, 3rd quantile: 106.5, max: 253647

saveRDS(ps.rare, here::here("Output/02 - Phyloseq Preprocessing Output/ps.rare.rds"))


