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

### Remove _sp.
as.data.frame(ps.All@tax_table) -> taxTable
gsub("_sp.", " Genus", taxTable$Species) -> taxTable.fix
as.data.frame(taxTable.fix) -> taxTable.fix

### Remove _
gsub("_", " ", taxTable.fix$taxTable.fix) -> taxTable.fix2
as.data.frame(taxTable.fix2) -> taxTable.fix2
taxTable$Taxa <- taxTable.fix2

tax_table(ps.All) <- as.matrix(taxTable)

## Convert phyloseq object to tree summarized experiment
tse <- mia::makeTreeSummarizedExperimentFromPhyloseq(ps.All)

## Set order for date and location
tse$Date_Long <- factor(tse$Date_Long, levels = c("August 2018 (Before)", "March 2019 (During)", "August 2019 (After)", "March 2020 (1 year later)"))


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
df_date_dunn$Species <- factor(df_date_dunn$Species, levels = c("Corallicola aquarius", "Cladocopium Genus", "Symbiodinium microadriaticum", "Other"))

## Plot
vol_plot_date <- df_date_dunn %>%
  ggplot(aes(x = LFC,
             y = -log10(adj.p),
             color = Species)) + 
  geom_point(size = 3.5, alpha = 0.8) + scale_color_manual(values = c("#FFB547", "#800000", "#FE9586", "darkgrey")) +
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



ggplot2::ggsave(here::here("plot_diffabund_date.png"), vol_plot_date,
                height = 250, width = 500, units = "mm",
                scale = 0.5, dpi = 1000)

