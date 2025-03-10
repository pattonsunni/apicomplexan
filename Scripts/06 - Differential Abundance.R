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


## A18 ====
subset_samples(ps.All, Date == "A18") -> ps.A18

## Convert phyloseq object to tree summarized experiment
tse <- mia::makeTreeSummarizedExperimentFromPhyloseq(ps.A18)

## Set order for location
tse$Location <- factor(tse$Location, levels = c("North", "East", "West"))

## Convert phyloseq object to tree summarized experiment
tse <- mia::makeTreeSummarizedExperimentFromPhyloseq(ps.A18)

## ANCOMBC2 formula
output_A18_ASV <- ancombc2(data = tse, assay_name = "counts", tax_level = "TaxonID", fix_formula = "Location", 
                           p_adj_method = "fdr", group = "Location", 
                           alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE, pairwise = TRUE,
                           mdfdr_control = list(fwer_ctrl_method = "fdr", B =100)) 


res_pair <- output_A18_ASV$res_pair
df_WE <- data.frame(c(res_pair[1], res_pair[4], res_pair[16], res_pair[19]), res_pair[22])
df_WE$Date <- "August 2018"
df_WE$Comp <- "West vs. East"
colnames(df_WE)[2] <- "LFC"
colnames(df_WE)[3] <- "adj.p"
colnames(df_WE)[4] <- "Diff"
colnames(df_WE)[5] <- "Pass_Sens"

df_NE <- data.frame(c(res_pair[1], res_pair[2], res_pair[14], res_pair[17]), res_pair[20])
df_NE$Date <- "August 2018"
df_NE$Comp <- "North vs. East"
colnames(df_NE)[2] <- "LFC"
colnames(df_NE)[3] <- "adj.p"
colnames(df_NE)[4] <- "Diff"
colnames(df_NE)[5] <- "Pass_Sens"

df_NW <- data.frame(c(res_pair[1], res_pair[3], res_pair[15], res_pair[18]), res_pair[21])
df_NW$Date <- "August 2018"
df_NW$Comp <- "North vs. West"
colnames(df_NW)[2] <- "LFC"
colnames(df_NW)[3] <- "adj.p"
colnames(df_NW)[4] <- "Diff"
colnames(df_NW)[5] <- "Pass_Sens"

## Merge dataframes
df_A18_pair <- dplyr::bind_rows(df_NW, df_NE, df_WE)

### Add another taxon column for plot coloring
df_A18_pair$Species <- paste0(df_A18_pair$taxon)
df_A18_pair$Species[df_A18_pair$adj.p >= 0.05] <- "Other"

subset(df_A18_pair, !(Diff == "TRUE" & Pass_Sens == "FALSE")) -> df_A18_pair

### Add another taxon column for plot shape
df_A18_pair$TaxNoASV <- paste0(df_A18_pair$Species)

gsub("ASV\\d+", "", df_A18_pair$TaxNoASV) -> df_A18_pair$TaxNoASV



vol_plot_date <- df_A18_pair %>%
  ggplot(aes(x = LFC,
             y = -log10(adj.p),
             color = TaxNoASV)) + 
  geom_point(size = 3.5, alpha = 0.8) +
  facet_wrap(~Date, labeller = labeller(Date_Long = label_wrap_gen(width = 15)))

vol_plot_date <- vol_plot_date + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed")

vol_plot_date + theme_bw(base_line_size = 1, base_rect_size = 1.5) +
  theme(axis.text = element_text(face = "bold", size = 14), 
        axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) +
  xlab("Log Fold Change") + ylab("-log10 Adjusted p value") + labs(color = 'Taxa') + theme(strip.text = element_text(face = "bold", size = 12)) 


## M19 ====
subset_samples(ps.All, Date == "M19") -> ps.M19

## Convert phyloseq object to tree summarized experiment
tse <- mia::makeTreeSummarizedExperimentFromPhyloseq(ps.M19)

## Set order for location
tse$Location <- factor(tse$Location, levels = c("North", "East", "West"))

## ANCOMBC2 formula
output_M19_ASV <- ancombc2(data = tse, assay_name = "counts", tax_level = "TaxonID", fix_formula = "Location", 
                           p_adj_method = "fdr", group = "Location", 
                           alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE, pairwise = TRUE,
                           mdfdr_control = list(fwer_ctrl_method = "fdr", B =100)) 


res_pair <- output_M19_ASV$res_pair
df_WE <- data.frame(c(res_pair[1], res_pair[4], res_pair[16], res_pair[19]), res_pair[22])
df_WE$Date <- "March 2019"
df_WE$Comp <- "West vs. East"
colnames(df_WE)[2] <- "LFC"
colnames(df_WE)[3] <- "adj.p"
colnames(df_WE)[4] <- "Diff"
colnames(df_WE)[5] <- "Pass_Sens"

df_NE <- data.frame(c(res_pair[1], res_pair[2], res_pair[14], res_pair[17]), res_pair[20])
df_NE$Date <- "March 2019"
df_NE$Comp <- "North vs. East"
colnames(df_NE)[2] <- "LFC"
colnames(df_NE)[3] <- "adj.p"
colnames(df_NE)[4] <- "Diff"
colnames(df_NE)[5] <- "Pass_Sens"

df_NW <- data.frame(c(res_pair[1], res_pair[3], res_pair[15], res_pair[18]), res_pair[21])
df_NW$Date <- "March 2019"
df_NW$Comp <- "North vs. West"
colnames(df_NW)[2] <- "LFC"
colnames(df_NW)[3] <- "adj.p"
colnames(df_NW)[4] <- "Diff"
colnames(df_NW)[5] <- "Pass_Sens"

## Merge dataframes
df_M19_pair <- dplyr::bind_rows(df_NW, df_NE, df_WE)

### Add another taxon column for plot coloring
df_M19_pair$Species <- paste0(df_M19_pair$taxon)
df_M19_pair$Species[df_M19_pair$adj.p >= 0.05] <- "Other"

subset(df_M19_pair, !(Diff == "TRUE" & Pass_Sens == "FALSE")) -> df_M19_pair

### Add another taxon column for plot shape
df_M19_pair$TaxNoASV <- paste0(df_M19_pair$Species)

gsub("ASV\\d+", "", df_M19_pair$TaxNoASV) -> df_M19_pair$TaxNoASV



vol_plot_date <- df_M19_pair %>%
  ggplot(aes(x = LFC,
             y = -log10(adj.p),
             color = TaxNoASV)) + 
  geom_point(size = 3.5, alpha = 0.8) + facet_grid(rows = vars(Date), cols = vars(Comp))
#facet_wrap(~Date, labeller = labeller(Date_Long = label_wrap_gen(width = 15)))

vol_plot_date <- vol_plot_date + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed")

vol_plot_date + theme_bw(base_line_size = 1, base_rect_size = 1.5) +
  theme(axis.text = element_text(face = "bold", size = 14), 
        axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) +
  xlab("Log Fold Change") + ylab("-log10 Adjusted p value") + labs(color = 'Taxa') + theme(strip.text = element_text(face = "bold", size = 12)) 


## A19 ====
subset_samples(ps.All, Date == "A19") -> ps.A19

## Convert phyloseq object to tree summarized experiment
tse <- mia::makeTreeSummarizedExperimentFromPhyloseq(ps.A19)
## Set order for location
tse$Location <- factor(tse$Location, levels = c("North", "East", "West"))

## ANCOMBC2 formula
output_A19_ASV <- ancombc2(data = tse, assay_name = "counts", tax_level = "TaxonID", fix_formula = "Location", 
                           p_adj_method = "fdr", group = "Location", 
                           alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE, pairwise = TRUE,
                           mdfdr_control = list(fwer_ctrl_method = "fdr", B =100)) 


res_pair <- output_A19_ASV$res_pair
df_WE <- data.frame(c(res_pair[1], res_pair[4], res_pair[16], res_pair[19]), res_pair[22])
df_WE$Date <- "August 2019"
df_WE$Comp <- "West vs. East"
colnames(df_WE)[2] <- "LFC"
colnames(df_WE)[3] <- "adj.p"
colnames(df_WE)[4] <- "Diff"
colnames(df_WE)[5] <- "Pass_Sens"

df_NE <- data.frame(c(res_pair[1], res_pair[2], res_pair[14], res_pair[17]), res_pair[20])
df_NE$Date <- "August 2019"
df_NE$Comp <- "North vs. East"
colnames(df_NE)[2] <- "LFC"
colnames(df_NE)[3] <- "adj.p"
colnames(df_NE)[4] <- "Diff"
colnames(df_NE)[5] <- "Pass_Sens"

df_NW <- data.frame(c(res_pair[1], res_pair[3], res_pair[15], res_pair[18]), res_pair[21])
df_NW$Date <- "August 2019"
df_NW$Comp <- "North vs. West"
colnames(df_NW)[2] <- "LFC"
colnames(df_NW)[3] <- "adj.p"
colnames(df_NW)[4] <- "Diff"
colnames(df_NW)[5] <- "Pass_Sens"

## Merge dataframes
df_A19_pair <- dplyr::bind_rows(df_NW, df_NE, df_WE)

### Add another taxon column for plot coloring
df_A19_pair$Species <- paste0(df_A19_pair$taxon)
df_A19_pair$Species[df_A19_pair$adj.p >= 0.05] <- "Other"

subset(df_A19_pair, !(Diff == "TRUE" & Pass_Sens == "FALSE")) -> df_A19_pair

### Add another taxon column for plot shape
df_A19_pair$TaxNoASV <- paste0(df_A19_pair$Species)

gsub("ASV\\d+", "", df_A19_pair$TaxNoASV) -> df_A19_pair$TaxNoASV



vol_plot_date <- df_A19_pair %>%
  ggplot(aes(x = LFC,
             y = -log10(adj.p),
             color = TaxNoASV)) + 
  geom_point(size = 3.5, alpha = 0.8) + facet_grid(cols = vars(Comp), rows = vars(Date))
#facet_wrap(~Date_Long, labeller = labeller(Date_Long = label_wrap_gen(width = 15)))

vol_plot_date <- vol_plot_date + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed")

vol_plot_date + theme_bw(base_line_size = 1, base_rect_size = 1.5) +
  theme(axis.text = element_text(face = "bold", size = 14), 
        axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) +
  xlab("Log Fold Change") + ylab("-log10 Adjusted p value") + labs(color = 'Taxa') + theme(strip.text = element_text(face = "bold", size = 12)) 



## M20 ====
subset_samples(ps.All, Date == "M20") -> ps.M20

## Convert phyloseq object to tree summarized experiment
tse <- mia::makeTreeSummarizedExperimentFromPhyloseq(ps.M20)

## Set order for location
tse$Location <- factor(tse$Location, levels = c("North", "East", "West"))

## ANCOMBC2 formula
output_M20_ASV <- ancombc2(data = tse, assay_name = "counts", tax_level = "TaxonID", fix_formula = "Location", 
                           p_adj_method = "fdr", group = "Location", 
                           alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE, pairwise = TRUE,
                           mdfdr_control = list(fwer_ctrl_method = "fdr", B =100)) 


res_pair <- output_M20_ASV$res_pair
df_WE <- data.frame(c(res_pair[1], res_pair[4], res_pair[16], res_pair[19]), res_pair[22])
df_WE$Date <- "March 2020"
df_WE$Comp <- "West vs. East"
colnames(df_WE)[2] <- "LFC"
colnames(df_WE)[3] <- "adj.p"
colnames(df_WE)[4] <- "Diff"
colnames(df_WE)[5] <- "Pass_Sens"

df_NE <- data.frame(c(res_pair[1], res_pair[2], res_pair[14], res_pair[17]), res_pair[20])
df_NE$Date <- "March 2020"
df_NE$Comp <- "North vs. East"
colnames(df_NE)[2] <- "LFC"
colnames(df_NE)[3] <- "adj.p"
colnames(df_NE)[4] <- "Diff"
colnames(df_NE)[5] <- "Pass_Sens"

df_NW <- data.frame(c(res_pair[1], res_pair[3], res_pair[15], res_pair[18]), res_pair[21])
df_NW$Date <- "March 2020"
df_NW$Comp <- "North vs. West"
colnames(df_NW)[2] <- "LFC"
colnames(df_NW)[3] <- "adj.p"
colnames(df_NW)[4] <- "Diff"
colnames(df_NW)[5] <- "Pass_Sens"

## Merge dataframes
df_M20_pair <- dplyr::bind_rows(df_NW, df_NE, df_WE)

### Add another taxon column for plot coloring
df_M20_pair$Species <- paste0(df_M20_pair$taxon)
df_M20_pair$Species[df_M20_pair$adj.p >= 0.05] <- "Other"

subset(df_M20_pair, !(Diff == "TRUE" & Pass_Sens == "FALSE")) -> df_M20_pair

### Add another taxon column for plot shape
df_M20_pair$TaxNoASV <- paste0(df_M20_pair$Species)

gsub("ASV\\d+", "", df_M20_pair$TaxNoASV) -> df_M20_pair$TaxNoASV



vol_plot_date <- df_M20_pair %>%
  ggplot(aes(x = LFC,
             y = -log10(adj.p),
             color = TaxNoASV)) + 
  geom_point(size = 3.5, alpha = 0.8) + facet_grid(cols = vars(Comp), rows = vars(Date))
#facet_wrap(~Date_Long, labeller = labeller(Date_Long = label_wrap_gen(width = 15)))

vol_plot_date <- vol_plot_date + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed")

vol_plot_date + theme_bw(base_line_size = 1, base_rect_size = 1.5) +
  theme(axis.text = element_text(face = "bold", size = 14), 
        axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) +
  xlab("Log Fold Change") + ylab("-log10 Adjusted p value") + labs(color = 'Taxa') + theme(strip.text = element_text(face = "bold", size = 12)) 


## Combine dfs ====
df_pairwise <- dplyr::bind_rows(df_M19_pair, df_A19_pair, df_M20_pair)

df_pairwise$Date <- factor(df_pairwise$Date, levels = c("March 2019", "August 2019", "March 2020"))

df_pairwise$TaxNoASV[df_pairwise$TaxNoASV == "Cladocopium  "] <- "Cladocopium"
df_pairwise$TaxNoASV[df_pairwise$TaxNoASV == "Cladocopium "] <- "Cladocopium"
df_pairwise$TaxNoASV[df_pairwise$TaxNoASV == "Symbiodinium "] <- "Symbiodinium"
df_pairwise$TaxNoASV[df_pairwise$TaxNoASV == "Corallicola aquarius "] <- "Corallicola aquarius"


df_pairwise$TaxNoASV <- factor(df_pairwise$TaxNoASV, levels = c("Corallicola aquarius", "Cladocopium", "Symbiodinium", "Other"))

colors <- c("#FFB547", "#800000", "#928EBE" , "lightgrey")


vol_plot_all <- df_pairwise %>%
  ggplot(aes(x = LFC,
             y = -log10(adj.p),
             color = TaxNoASV)) + 
  geom_point(size = 3.5, alpha = 0.8) + facet_grid(cols = vars(Comp), rows = vars(Date)) + scale_color_manual(values = colors)

vol_plot_all <- vol_plot_all + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed")

vol_plot_all <- vol_plot_all + theme_bw(base_line_size = 1, base_rect_size = 1.5) +
  theme(axis.text = element_text(face = "bold", size = 14), 
        axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) +
  xlab("Log Fold Change") + ylab("-log10 Adjusted p value") + labs(color = 'Taxa') + theme(strip.text = element_text(face = "bold", size = 12)) 


write.csv(df_pairwise, here::here("Output/06 - Differential Abundance Output - DiffAbund_ASV_compare_shore.csv"))

ggplot2::ggsave(here::here("Output/06 - Differential Abundance Output/plot_DA_shore.png"), vol_plot_all,
                height = 400, width = 600, units = "mm",
                scale = 0.5, dpi = 1000)





## Below is old (doing dunnett procedure)

# subset
subset_samples(ps.All, Date == "A19") -> ps.A19

## Convert phyloseq object to tree summarized experiment
tse <- mia::makeTreeSummarizedExperimentFromPhyloseq(ps.A19)

## Set order for date and location
tse$Date_Long <- factor(tse$Date_Long, levels = c("August 2018 (Before)", "March 2019 (During)", "August 2019 (After)", "March 2020 (1 year later)"))

## ANCOMBC2 formula
### CoralID as random effect
#output_east_ASV<- ancombc2(data = tse, assay_name = "counts", tax_level = "Taxa", fix_formula = "Date_Long", 
 #                             rand_formula = "(1 | CoralID)", p_adj_method = "fdr", group = "Date_Long", 
  #                            alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE, dunnet = TRUE,
   #                           mdfdr_control = list(fwer_ctrl_method = "fdr", B =100)) 

output_A19_ASV <- ancombc2(data = tse, assay_name = "counts", tax_level = "TaxonID", fix_formula = "Location", 
                                                        p_adj_method = "fdr", group = "Location", 
                                                        alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE, pairwise = TRUE,
                                                        mdfdr_control = list(fwer_ctrl_method = "fdr", B =100)) 
                           

## Save output 
res_dunn <- output_east_ASV$res_dunn

### Save specific output as dataframe 
#### M19
df_M19_dunn <- data.frame(c(res_dunn[1], res_dunn[2], res_dunn[14], res_dunn[17]), res_dunn[20])
df_M19_dunn$Date_Long <- "March 2019 (During)"
colnames(df_M19_dunn)[2] <- "LFC"
colnames(df_M19_dunn)[3] <- "adj.p"
colnames(df_M19_dunn)[4] <- "Diff"
colnames(df_M19_dunn)[5] <- "Pass_Sens"

#### A19
df_A19_dunn <- data.frame(c(res_dunn[1], res_dunn[3], res_dunn[15], res_dunn[18], res_dunn[21]))
df_A19_dunn$Date_Long <- "August 2019 (After)"
colnames(df_A19_dunn)[2] <- "LFC"
colnames(df_A19_dunn)[3] <- "adj.p"
colnames(df_A19_dunn)[4] <- "Diff"
colnames(df_A19_dunn)[5] <- "Pass_Sens"


#### M20
df_M20_dunn <- data.frame(c(res_dunn[1], res_dunn[4], res_dunn[16], res_dunn[19], res_dunn[22]))
df_M20_dunn$Date_Long <- "March 2020 (1 year later)"
colnames(df_M20_dunn)[2] <- "LFC"
colnames(df_M20_dunn)[3] <- "adj.p"
colnames(df_M20_dunn)[4] <- "Diff"
colnames(df_M20_dunn)[5] <- "Pass_Sens"


## Merge dataframes
df_date_dunn <- dplyr::bind_rows(df_M19_dunn, df_A19_dunn, df_M20_dunn)

### Add another taxon column for plot coloring
df_date_dunn$Species <- paste0(df_date_dunn$taxon)
df_date_dunn$Species[df_date_dunn$adj.p >= 0.05] <- "Other"

subset(df_date_dunn, !(Diff == "TRUE" & Pass_Sens == "FALSE")) -> df_date_dunn

### Add another taxon column for plot shape
df_date_dunn$TaxNoASV <- paste0(df_date_dunn$Species)

gsub("ASV\\d+", "", df_date_dunn$TaxNoASV) -> df_date_dunn$TaxNoASV


write.csv(df_date_dunn, here::here("ancombc2_output.csv"))

# Set order for plotting
df_date_dunn$Date_Long <- factor(df_date_dunn$Date_Long, levels = c("March 2019 (During)", "August 2019 (After)",
                                                               "March 2020 (1 year later)"))
#df_date_dunn$Species <- factor(df_date_dunn$Species, levels = c("Corallicola aquarius", "Cladocopium Genus", "Symbiodinium microadriaticum", "Other"))

## Plot
vol_plot_date <- df_date_dunn %>%
  ggplot(aes(x = LFC,
             y = -log10(adj.p),
             color = TaxNoASV)) + 
  geom_point(size = 3.5, alpha = 0.8) + #aes(shape = factor(TaxNoASV)), 
#+ scale_color_manual(values = c("#FFB547", "#800000", "#FE9586", "darkgrey")) +
facet_wrap(~Date_Long, labeller = labeller(Date_Long = label_wrap_gen(width = 15)))

vol_plot_date <- vol_plot_date + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed")
  # +
#  scale_x_continuous(breaks = c(seq(-6, 6, 2)), # Modify x-axis tick intervals    
 #                    limits = c(-7, 7)) #+(0, 5)

vol_plot_date <- vol_plot_date + theme_bw(base_line_size = 1, base_rect_size = 1.5) +
  theme(axis.text = element_text(face = "bold", size = 14), 
        axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) +
  xlab("Log Fold Change") + ylab("-log10 Adjusted p value") + labs(color = 'Taxa') + theme(strip.text = element_text(face = "bold", size = 12)) 



ggplot2::ggsave(here::here("plot_diffabund_date.png"), vol_plot_date,
                height = 250, width = 500, units = "mm",
                scale = 0.5, dpi = 1000)

## diff abund not at ASV level
## ANCOMBC2 formula
### CoralID as random effect
output_time <- ancombc2(data = tse, assay_name = "counts", tax_level = "Taxa", fix_formula = "Date_Long + Location", 
                        rand_formula = "(1 | CoralID)", p_adj_method = "fdr", group = "Date_Long", 
                        alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE, dunnet = TRUE,
                        mdfdr_control = list(fwer_ctrl_method = "fdr", B =100)) 


## Save output 
res_dunn <- output_time$res_dunn

### Save specific output as dataframe 
#### M19
df_M19_dunn <- data.frame(c(res_dunn[1], res_dunn[2], res_dunn[14], res_dunn[17]), res_dunn[20])
df_M19_dunn$Date_Long <- "March 2019 (During)"
colnames(df_M19_dunn)[2] <- "LFC"
colnames(df_M19_dunn)[3] <- "adj.p"
colnames(df_M19_dunn)[4] <- "Diff"
colnames(df_M19_dunn)[5] <- "Pass_Sens"

#### A19
df_A19_dunn <- data.frame(c(res_dunn[1], res_dunn[3], res_dunn[15], res_dunn[18], res_dunn[21]))
df_A19_dunn$Date_Long <- "August 2019 (After)"
colnames(df_A19_dunn)[2] <- "LFC"
colnames(df_A19_dunn)[3] <- "adj.p"
colnames(df_A19_dunn)[4] <- "Diff"
colnames(df_A19_dunn)[5] <- "Pass_Sens"


#### M20
df_M20_dunn <- data.frame(c(res_dunn[1], res_dunn[4], res_dunn[16], res_dunn[19], res_dunn[22]))
df_M20_dunn$Date_Long <- "March 2020 (1 year later)"
colnames(df_M20_dunn)[2] <- "LFC"
colnames(df_M20_dunn)[3] <- "adj.p"
colnames(df_M20_dunn)[4] <- "Diff"
colnames(df_M20_dunn)[5] <- "Pass_Sens"


## Merge dataframes
df_date_dunn <- dplyr::bind_rows(df_M19_dunn, df_A19_dunn, df_M20_dunn)

### Add another taxon column for plot coloring
df_date_dunn$Species <- paste0(df_date_dunn$taxon)
df_date_dunn$Species[df_date_dunn$adj.p >= 0.05] <- "Other"




subset(df_date_dunn, !(Diff == "TRUE" & Pass_Sens == "FALSE")) -> df_date_dunn

write.csv(df_date_dunn, here::here("ancombc2_output.csv"))

# Set order for plotting
df_date_dunn$Date_Long <- factor(df_date_dunn$Date_Long, levels = c("March 2019 (During)", "August 2019 (After)",
                                                                    "March 2020 (1 year later)"))
#df_date_dunn$Species <- factor(df_date_dunn$Species, levels = c("Corallicola aquarius", "Cladocopium Genus", "Symbiodinium microadriaticum", "Other"))


## Plot
vol_plot_date <- df_date_dunn %>%
  ggplot(aes(x = LFC,
             y = -log10(adj.p),
             color = Species)) + 
  geom_point(size = 3.5, alpha = 0.8) + #+ scale_color_manual(values = c("#FFB547", "#800000", "#FE9586", "darkgrey")) +
  facet_wrap(~Date_Long, labeller = labeller(Date_Long = label_wrap_gen(width = 15)))

vol_plot_date <- vol_plot_date + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed")# +
# +
#  scale_x_continuous(breaks = c(seq(-6, 6, 2)), # Modify x-axis tick intervals    
#                    limits = c(-7, 7)) #+(0, 5)

vol_plot_date <- vol_plot_date + theme_bw(base_line_size = 1, base_rect_size = 1.5) +
  theme(axis.text = element_text(face = "bold", size = 14), 
        axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) +
  xlab("Log Fold Change") + ylab("-log10 Adjusted p value") + labs(color = 'Taxa') + theme(strip.text = element_text(face = "bold", size = 12)) 
