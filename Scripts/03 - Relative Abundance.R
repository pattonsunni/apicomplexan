# Author: Sunni Patton
# Last edited: 11/08/2024
# Title: Relative abundance
# Overview: Assessing microeukaryote relative abundance

## Set seed ====
set.seed(123)

## Load libraries ====
library(ggplot2)
library(phyloseq)

## Load data ====
readRDS(here::here("Output/02 - Phyloseq Preprocessing Output/ps.rare.rds")) -> ps.rare

## Transform to relative abundance ====
ps.rare.trans <- transform_sample_counts(ps.rare, function(OTU) OTU/sum(OTU))

## Prune most abundant taxa ====
top100 <- names(sort(taxa_sums(ps.rare.trans), decreasing = TRUE))[1:100] 
ps.rare.top100 <- prune_taxa(top100, ps.rare.trans)

# Change coralID to factor
ps.rare.top100@sam_data$CoralID <- as.factor(ps.rare@sam_data$CoralID)

## Plot relative abundance for all non-metazoan micro-eukaryotes ====
tax_glom(ps.rare.top100, taxrank = "Species", NArm = FALSE) -> ps.rare.top100

plot.relAbund <- plot_bar(subset_samples(ps.rare.top100), x="CoralID", fill="Species") + 
  facet_grid(vars(Date_Long), vars(LTER_Site), scales = "free_x") 
plot.relAbund <- plot.relAbund + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 11.5), 
        axis.text.y = element_text(face = "bold", size = 11.5), title = element_text(face = "bold")) + xlab("CoralID") + labs(title = "Relative Abundance - Non-metazoan micro-eukaryotes")

## Subset only alveolates ====
subset_taxa(ps.rare, Division == "Alveolata") -> ps.rare.alv

## Transform to relative abundance ====
ps.rare.trans <- transform_sample_counts(ps.rare.alv, function(OTU) OTU/sum(OTU))

## Prune most abundant taxa ====
top100 <- names(sort(taxa_sums(ps.rare.trans), decreasing = TRUE))[1:100] 
ps.rare.top100 <- prune_taxa(top100, ps.rare.trans)

# Change coralID to factor
ps.rare.top100@sam_data$CoralID <- as.factor(ps.rare@sam_data$CoralID)

## Plot relative abundance for alveolates ====
tax_glom(ps.rare.top100, taxrank = "Species", NArm = FALSE) -> ps.rare.top100

plot.relAbund <- plot_bar(subset_samples(ps.rare.top100), x="CoralID", fill="Species") + 
  facet_grid(vars(Date_Long), vars(LTER_Site), scales = "free_x") 
plot.relAbund <- plot.relAbund + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 11.5), 
        axis.text.y = element_text(face = "bold", size = 11.5), title = element_text(face = "bold")) + xlab("CoralID") + labs(title = "Relative Abundance - Alveolates")

