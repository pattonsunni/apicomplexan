# Author: Sunni Patton
# Last edited: 11/22/24
# Title: Differential abundance
# Overview: Identifying differentially abundant taxa

## Set seed ====
set.seed(123)

## Load packages ====
library(dplyr)
library(ggplot2)
library(tidyr)
library(ANCOMBC)
library(phyloseq)
library(mia)
library(lme4)
library(lmerTest)

## Load data ====
readRDS(here::here("Output/02 - Phyloseq Preprocessing Output/ps.All.rds")) -> ps.All

## Make sure site is factor
ps.All@sam_data$LTER_Site <- as.factor(ps.All@sam_data$LTER_Site)

## Add column for long date
ps.All@sam_data$Date_Long <- paste0(ps.All@sam_data$Date)
### Change to long format
ps.All@sam_data$Date_Long[ps.All@sam_data$Date_Long == "A18"] <- "August 2018 (Before)"
ps.All@sam_data$Date_Long[ps.All@sam_data$Date_Long == "M19"] <- "March 2019 (During)"
ps.All@sam_data$Date_Long[ps.All@sam_data$Date_Long == "A19"] <- "August 2019 (After)"
ps.All@sam_data$Date_Long[ps.All@sam_data$Date_Long == "M20"] <- "March 2020 (1 year later)"

## Subset by shore, then look at date ====
ps.north <- subset_samples(ps.All, Location == "North")

## Convert phyloseq object to tree summarized experiment
tse <- mia::makeTreeSummarizedExperimentFromPhyloseq(ps.north)

## Set order for date and location
tse$Date_Long <- factor(tse$Date_Long, levels = c("August 2018 (Before)", "March 2019 (During)", "August 2019 (After)", "March 2020 (1 year later)"))
tse$Location <- factor(tse$Location, levels = c("North", "East", "West"))


## ANCOMBC2 formula
### CoralID as random effect
output_time <- ancombc2(data = tse, assay_name = "counts", tax_level = "Species", fix_formula = "Date_Long", 
                       rand_formula = "(1 | CoralID)", p_adj_method = "fdr", group = "Date_Long", 
                       alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE) 

output_location <- ancombc2(data = tse, assay_name = "counts", tax_level = "Species", fix_formula = "Location", 
                        rand_formula = "(1 | CoralID)", p_adj_method = "fdr", group = "Location", 
                        alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE) 


## Save output 
res_prim <- output_time$res

res_prim_loc <- output_location$res

### Save specific output as dataframe 
#### M19
df_M19 <- data.frame(c(res_prim[1], res_prim[3], res_prim[19], res_prim[23]))
df_M19$Date_Long <- "March 2019 (During)"
colnames(df_M19)[2] <- "LFC"
colnames(df_M19)[3] <- "adj.p"
colnames(df_M19)[4] <- "Diff"

#### A19
df_A19 <- data.frame(c(res_prim[1], res_prim[4], res_prim[20], res_prim[24]))
df_A19$Date_Long <- "August 2019 (After)"
colnames(df_A19)[2] <- "LFC"
colnames(df_A19)[3] <- "adj.p"
colnames(df_A19)[4] <- "Diff"

#### M20
df_M20 <- data.frame(c(res_prim[1], res_prim[5], res_prim[21], res_prim[25]))
df_M20$Date_Long <- "March 2020 (1 year later)"
colnames(df_M20)[2] <- "LFC"
colnames(df_M20)[3] <- "adj.p"
colnames(df_M20)[4] <- "Diff"

## Merge dataframes
df_date <- dplyr::bind_rows(df_M19, df_A19, df_M20)

## Save specific output as dataframe
df_east <- data.frame(c(res_prim_loc[1], res_prim_loc[3], res_prim_loc[15], res_prim_loc[18]))
df_east$Location <- "East"
colnames(df_east)[2] <- "LFC"
colnames(df_east)[3] <- "adj.p"
colnames(df_east)[4] <- "Diff"

df_west <- data.frame(c(res_prim_loc[1], res_prim_loc[4], res_prim_loc[16], res_prim_loc[19]))
df_west$Location <- "West"
colnames(df_west)[2] <- "LFC"
colnames(df_west)[3] <- "adj.p"
colnames(df_west)[4] <- "Diff"

df_loc <- bind_rows(df_east, df_west)

### Add another taxon column for plot coloring
df_date$Species <- paste0(df_date$taxon)
df_date$Species[df_date$adj.p >= 0.05] <- "Other"

df_loc$Species <- paste0(df_loc$taxon)
df_loc$Species[df_loc$adj.p >= 0.05] <- "Other"

# Set order for plotting
df_date$Date_Long <- factor(df_date$Date_Long, levels = c("March 2019 (During)", "August 2019 (After)",
                                                          "March 2020 (1 year later)"))

### Add another taxon column for plot coloring
df_date$Species2 <- paste0(df_date$Species)
df_date$Species2[df_date$adj.p >= 0.05] <- "Other"
df_date$Species2[df_date$Species2 == "Species:Cladocopium_sp."] <- "Cladocopium" 
df_date$Species2[df_date$Species2 == "Species:Corallicola_aquarius"] <- "Corallicola aquarius" 
df_date$Species2[df_date$Species2 == "Species:Corallicola_sp."] <- "Corallicola sp." 
df_date$Species2[df_date$Species2 == "Species:Licnophora_macfarlandi"] <- "Licnophora macfarlandi"
df_date$Species2[df_date$Species2 == "Species:Symbiodinium_microadriaticum"] <- "Symbiodinium microadriaticum"
df_date$Species2[df_date$Species2 == "Species:Symbiodinium_sp."] <- "Symbiodinium sp."
df_date$Species2[df_date$Species2 == "Genus:Symbiodinium"] <- "Symbiodinium"

df_date$Species2 <- factor(df_date$Species2, levels = c("Cladocopium", "Corallicola aquarius", 
                                                      "Corallicola sp.", "Licnophora macfarlandi", 
                                                      "Symbiodinium microadriaticum", "Symbiodinium sp.",
                                                      "Symbiodinium", "Other"))

## Plot
vol_plot_date <- df_date %>%
  ggplot(aes(x = LFC,
             y = -log10(adj.p),
             color = Species)) + 
  geom_point(size = 3.5, alpha = 0.8) + scale_color_manual(values = c("#800000", "#FFB547",
                                                                      "#ADB17D","#155F83", 
                                                                      "#C16622", "#928EBE", "#D25F9E", "darkgrey")) +facet_wrap(~Date_Long, labeller = labeller(Date_Long = label_wrap_gen(width = 15)))

#"#70AF94", "#6F99AD","#928EBE","#D25F9E", "#FEC655", "#B24745", "#00A1D5", "darkgrey"

vol_plot_date <- vol_plot_date + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  scale_x_continuous(breaks = c(seq(-6, 6, 2)), # Modify x-axis tick intervals    
                     limits = c(-7, 7)) #+(0, 5)

vol_plot_date <- vol_plot_date + theme_bw(base_line_size = 1, base_rect_size = 1.5) +
  theme(axis.text = element_text(face = "bold", size = 14), 
        axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) +
  xlab("Log Fold Change") + ylab("-log10 Adjusted p value") + labs(color = 'Species') + theme(strip.text = element_text(face = "bold", size = 12)) 



### Add another taxon column for plot coloring
df_loc$Species2 <- paste0(df_loc$Species)
df_loc$Species2[df_loc$adj.p >= 0.05] <- "Other"
df_loc$Species2[df_loc$Species2 == "Species:Cladocopium_sp."] <- "Cladocopium" 
df_loc$Species2[df_loc$Species2 == "Species:Corallicola_aquarius"] <- "Corallicola aquarius" 
df_loc$Species2[df_loc$Species2 == "Family:Corallicolidae"] <- "Corallicolidae" 
df_loc$Species2[df_loc$Species2 == "Genus:Ceramium"] <- "Ceramium"
df_loc$Species2[df_loc$Species2 == "Species:Symbiodinium_sp."] <- "Symbiodinium sp."

df_loc$Species2 <- factor(df_loc$Species2, levels = c("Ceramium", "Cladocopium", "Corallicola aquarius", 
                                                        "Corallicolidae", "Symbiodinium sp.", "Other"))


vol_plot_loc <- df_loc %>%
  ggplot(aes(x = LFC,
             y = -log10(adj.p),
             color = Species2)) + 
  geom_point(size = 3.5, alpha = 0.8) + scale_color_manual(values = c("#075149", "#800000", "#FFB547", "#FD7446", "#928EBE", "darkgrey")) +facet_wrap(~Location)

vol_plot_loc <- vol_plot_loc + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  scale_x_continuous(breaks = c(seq(-6, 6, 2)), # Modify x-axis tick intervals    
                     limits = c(-7, 7)) #+(0, 5)

vol_plot_loc <- vol_plot_loc + theme_bw(base_line_size = 1, base_rect_size = 1.5) +
  theme(axis.text = element_text(face = "bold", size = 14), 
        axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) +
  xlab("Log Fold Change") + ylab("-log10 Adjusted p value") + labs(color = 'Species') + theme(strip.text = element_text(face = "bold", size = 12)) 



ggplot2::ggsave(here::here("Output/06 - Differential Abundance Output/plot_diffabund_date.png"), vol_plot_date,
                height = 250, width = 500, units = "mm",
                scale = 0.5, dpi = 1000)

ggplot2::ggsave(here::here("Output/06 - Differential Abundance Output/plot_diffabund_loc.png"), vol_plot_loc,
                height = 250, width = 400, units = "mm",
                scale = 0.5, dpi = 1000)
