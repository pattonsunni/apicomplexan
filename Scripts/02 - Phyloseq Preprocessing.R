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
library(microViz)
library(ggplot2)

## Load data ====
read.csv(here::here("Output/01 - Metadata Output/metadata.csv")) -> metadata
readRDS(here::here("Output/00 - Read Preprocessing Output/seq_table_nometa.rds")) -> seq_table_nometa
readRDS(here::here("Output/00 - Read Preprocessing Output/taxa_nometa.rds")) -> taxa_nometa

read.csv(here::here("Quant.csv")) -> quant.df
merge(metadata, quant.df) -> metadata

## Edit sample data ====
### Make sample numbers into sample names to work with phyloseq
x <- metadata$SampleID
rownames(metadata) <- x

## Create phyloseq object ====
ps.All <- phyloseq(otu_table(seq_table_nometa, taxa_are_rows=FALSE), sample_data(metadata), tax_table(taxa_nometa))

## Remove contaminants ====
sample_data(ps.All)$Sample_or_Control <- "Sample"
sample_data(ps.All)$Sample_or_Control[sample_data(ps.All)$LTER_Site == "NC"] <- "Negative Control"
sample_data(ps.All)$is.neg <- sample_data(ps.All)$Sample_or_Control == "Negative Control"

### Plot sample library size 
df <- as.data.frame(sample_data(ps.All)) 
df$LibrarySize <- sample_sums(ps.All)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
plot_decontam <- ggplot(data=df, aes(x = Index, y = LibrarySize, color = is.neg)) + geom_point()

### Identify contaminants
contamdf_combined <- isContaminant(ps.All, neg="is.neg", conc="quant_reading", method="combined", threshold=0.4)
table(contamdf_combined$contaminant) #identified 4 as contaminants 
head(which(contamdf_combined$contaminant)) 

ps.All <- prune_taxa(!contamdf_combined$contaminant, ps.All) # Contains 1212 ASVs after contaminants removed

subset_taxa(ps.All, Domain != "NA") -> ps.All

## Fix taxonomy table ====
ps.All <- phyloseq_validate(ps.All)

tax_fix(
  ps.All,
  min_length = 4,
  unknowns = NA,
  suffix_rank = "classified",
  sep = " ",
  anon_unique = TRUE,
  verbose = TRUE
) -> ps.All

## Remove NCs
subset_samples(ps.All, Location != "NC") -> ps.All # 85 samples
ps.All <- prune_taxa(taxa_sums(ps.All@otu_table) > 0, ps.All) # 1209 total taxa

## Make rarefaction curve ====
as.matrix(as.data.frame(ps.All@otu_table)) -> data
sort(sample_sums(ps.All))


# Make rarecurve to see how much diversity is captured if we rarefy to lowest sample depth (1084)
vegan::rarecurve(data, step = 100, label = FALSE, abline(v = 1084), col = "blue") # Essentially no diversity is captured
# Make rarecurve to see how much diversity is captured if we rarefy to 14205
vegan::rarecurve(data, step = 100, label = FALSE, abline(v = 14205, col = "red"), col = "black") # Most diversity is captured

## If we choose to rarefy to 14205, we automatically lose 2 samples and we still won't capture all the diversity; but will get a majority

## Remove samples with ~1000 reads ====
subset_samples(ps.All, SampleID != "1_BAK_ACR_57_M19") -> ps.All
subset_samples(ps.All, SampleID != "1_BAK_ACR_57_A18") -> ps.All

# Save phyloseq object
ps.All <- prune_taxa(taxa_sums(ps.All@otu_table) > 0, ps.All) # 1205 total taxa
saveRDS(ps.All, here::here("ps.All.rds"))  
sum(sample_sums(ps.All)) # 10319623 total reads

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
ps.rare <- prune_taxa(taxa_sums(ps.rare@otu_table) > 0, ps.rare) # Removed 187 taxa (1018 remaining)

summary(taxa_sums(ps.rare)) # minimum times a taxa appears: 1, 1st quartile: 4, median: 19.0, mean: 1158.2     3rd quantile: 100.8 , max: 253634.0 

saveRDS(ps.rare, here::here("ps.rare.rds"))


