# Author: Sunni Patton
# Last edited: 11/08/2024
# Title: Relative abundance
# Overview: Assessing microeukaryote relative abundance

## Set seed ====
set.seed(123)

## Load libraries ====
library(ggplot2)
library(phyloseq)
library(speedyseq)
library(microViz)
library(ggpubr)

## Load data ====
readRDS(here::here("Output/02 - Phyloseq Preprocessing Output/ps.rare.rds")) -> ps.rare

## Validate phyloseq object and fix if necessary ====
ps.rare <- phyloseq_validate(ps.rare)

tax_fix(
  ps.rare,
  min_length = 4,
  unknowns = NA,
  suffix_rank = "classified",
  sep = " ",
  anon_unique = TRUE,
  verbose = TRUE
) -> ps.rare


## Transform to relative abundance ====
ps.rare.trans <- transform_sample_counts(ps.rare, function(OTU) OTU/sum(OTU))

## Prune most abundant taxa ====
top100 <- names(sort(taxa_sums(ps.rare.trans), decreasing = TRUE))[1:100] 
ps.rare.top100 <- prune_taxa(top100, ps.rare.trans)

# Change coralID to factor
ps.rare.top100@sam_data$CoralID <- as.factor(ps.rare@sam_data$CoralID)

## Plot relative abundance for all non-metazoan micro-eukaryotes ====
tax_glom(ps.rare.top100, taxrank = "Species", NArm = FALSE) -> ps.rare.top100

colors <- c("#4BACC6", "#EA6312", "#800000", "#FFB547", "#ADB17D", "#7F5F52", "#D092A7",
            "#9B6BF2", "#E33D6F", "#8EC0C1", "#4472C4", "#928EBE", "lightgrey")

as.data.frame(ps.rare.top100@tax_table) -> x
x$Species[x$Species == "Cladocopium_sp."] <- "Cladocopium sp." 
x$Species[x$Species == "Bacillariophyceae Class"] <- "Bacillariophyceae Family" 
x$Species[x$Species == "Cladocopium Genus"] <- "Cladocopium sp." 
x$Species[x$Species == "Symbiodinium_sp."] <- "Symbiodinium sp." 
x$Species[x$Species == "Corallicola_aquarius"] <- "Corallicola aquarius" 
x$Species[x$Species == "Corallicola_sp."] <- "Corallicola sp." 
x$Species[x$Species == "Cylindrotheca_closterium"] <- "Cylindrotheca closterium" 
x$Species[x$Species == "Hydrolithon_onkodes"] <- "Hydrolithon onkodes" 
x$Species[x$Species == "Polysiphonia_sp."] <- "Polysiphonia sp." 
x$Species[x$Species == "Rhodomelaceae_X_sp."] <- "Rhodomelaceae Family" 
x$Species[x$Species == "Bacillaiophyceae Class"] <- "Bacillaiophyceae Family" 


tax_table(ps.rare.top100) <- as.matrix(x)

ps.rare.top100@sam_data$Date_Long <- factor(ps.rare.top100@sam_data$Date_Long, levels = c("August 2018", "March 2019", "August 2019","March 2020"))


plot.relAbund <- plot_bar(ps.rare.top100, x="CoralID", fill="Species") + 
  facet_grid(vars(Date_Long), vars(LTER_Site), scales = "free_x") + scale_fill_manual(values = colors)
plot.relAbund <- plot.relAbund + theme_bw(base_line_size = 1, base_rect_size = 1.5) + 
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 11.5), 
        axis.text.y = element_text(face = "bold", size = 11.5), title = element_text(face = "bold")) + xlab("Coral ID") + ylab("Relative Abundance (%)") +
  labs(title = "Non-metazoan micro-eukaryotes") + theme(strip.text = element_text(face = "bold", size = 11)) 

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

as.data.frame(ps.rare.top100@tax_table) -> x
x$Species[x$Species == "Cladocopium_sp."] <- "Cladocopium sp." 
x$Species[x$Species == "Cladocopium Genus"] <- "Cladocopium sp." 
x$Species[x$Species == "Corallicola_aquarius"] <- "Corallicola aquarius" 
x$Species[x$Species == "Corallicola_sp."] <- "Corallicola sp." 
x$Species[x$Species == "Symbiodinium_sp."] <- "Symbiodinium sp." 
x$Species[x$Species == "Protocruzia_contrax"] <- "Protocruzia contrax" 
x$Species[x$Species == "Coolia_canariensis"] <- "Coolia canariensis" 

tax_table(ps.rare.top100) <- as.matrix(x)

colors <- c("#800000", "#A1D68B", "#FFB547", "#ADB17D", "#7F5F52", "#B31166", "#C6C1F0", "#4472C4", "#928EBE", "lightgrey")

ps.rare.top100@sam_data$Date_Long <- factor(ps.rare.top100@sam_data$Date_Long, levels = c("August 2018", "March 2019", "August 2019","March 2020"))


plot.relAbund_alv <- plot_bar(subset_samples(ps.rare.top100), x="CoralID", fill="Species") + 
  facet_grid(vars(Date_Long), vars(LTER_Site), scales = "free_x") + scale_fill_manual(values = colors)
plot.relAbund_alv <- plot.relAbund_alv + theme_bw(base_line_size = 1, base_rect_size = 1.5) + 
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 11.5), 
        axis.text.y = element_text(face = "bold", size = 11.5), title = element_text(face = "bold")) + 
  xlab("Coral ID") + ylab("Relative Abundance (%)") + labs(title = "Alveolates") + theme(strip.text = element_text(face = "bold", size = 11))

plots <- ggarrange(plot.relAbund, plot.relAbund_alv, labels = c("A", "B"), ncol = 1, nrow = 2)

ggplot2::ggsave(here::here("Output/03 - Relative Abundance Output/plot_relabund_new.png"), plots,
                height = 600, width = 500, units = "mm",
                scale = 0.5, dpi = 1000)
