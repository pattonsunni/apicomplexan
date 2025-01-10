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
library(dplyr)
library(rstatix)
library(broom)

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

### Remove _sp.
as.data.frame(ps.rare@tax_table) -> taxTable
gsub("_sp.", " Genus", taxTable$Species) -> taxTable.fix
as.data.frame(taxTable.fix) -> taxTable.fix

### Remove _
gsub("_", " ", taxTable.fix$taxTable.fix) -> taxTable.fix2
as.data.frame(taxTable.fix2) -> taxTable.fix2
taxTable$Taxa <- taxTable.fix2

tax_table(ps.rare) <- as.matrix(taxTable)

saveRDS(ps.rare, here::here("ps.rare.renamed.rds"))

## Transform to relative abundance ====
ps.rare.trans <- transform_sample_counts(ps.rare, function(OTU) OTU/sum(OTU))

## Prune most abundant taxa ====
top100 <- names(sort(taxa_sums(ps.rare.trans), decreasing = TRUE))[1:100] 
ps.rare.top100 <- prune_taxa(top100, ps.rare.trans)

# Change coralID to factor
ps.rare.top100@sam_data$CoralID <- as.factor(ps.rare@sam_data$CoralID)

## Plot relative abundance for all non-metazoan micro-eukaryotes ====
tax_glom(ps.rare.top100, taxrank = "Taxa", NArm = FALSE) -> ps.rare.top100

colors <- c("#4BACC6", "#EA6312", "#800000", "#FFB547", "#ADB17D", "#7F5F52", "#D092A7",
            "#9B6BF2", "#E33D6F", "#8EC0C1", "#4472C4", "#928EBE", "#24693D")

as.data.frame(ps.rare.top100@tax_table) -> x
x$Taxa[x$Taxa == "Rhodomelaceae X Genus"] <- "Rhodomelaceae Family" 

tax_table(ps.rare.top100) <- as.matrix(x)

ps.rare.top100@sam_data$Date_Long <- factor(ps.rare.top100@sam_data$Date_Long, levels = c("August 2018", "March 2019", "August 2019","March 2020"))


## Plot by site and time
plot.relAbund <- plot_bar(ps.rare.top100, x="CoralID", fill="Taxa") + 
  facet_grid(vars(Date_Long), vars(LTER_Site), scales = "free_x") + scale_fill_manual(values = colors)
plot.relAbund <- plot.relAbund + theme_bw(base_line_size = 1, base_rect_size = 1.5) + 
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 11.5), 
        axis.text.y = element_text(face = "bold", size = 11.5), title = element_text(face = "bold")) + xlab("Coral ID") + ylab("Relative Abundance (%)") +
  labs(title = "Non-metazoan microeukaryotes") + theme(strip.text = element_text(face = "bold", size = 11)) + labs(fill = "Taxa") 


## Plot by location and time
ps.rare.top100@sam_data$Location <- factor(ps.rare.top100@sam_data$Location, levels = c("North", "East", "West"))

plot.relAbund_loc <- plot_bar(ps.rare.top100, x="CoralID", fill="Taxa") + 
  facet_grid(vars(Date_Long), vars(Location), scales = "free_x") + scale_fill_manual(values = colors)
plot.relAbund_loc <- plot.relAbund_loc + theme_bw(base_line_size = 1, base_rect_size = 1.5) + 
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 11.5), 
        axis.text.y = element_text(face = "bold", size = 11.5), title = element_text(face = "bold")) + xlab("Coral ID") + ylab("Relative Abundance (%)") +
  labs(title = "Non-metazoan microeukaryotes") + theme(strip.text = element_text(face = "bold", size = 11)) + labs(fill = "Taxa") 



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
tax_glom(ps.rare.top100, taxrank = "Taxa", NArm = FALSE) -> ps.rare.top100

colors <- c("#800000","#A1D68B", "#FFB547", "#ADB17D", "#7F5F52", "#B31166", "#24693D", "#F3E79A", "#4472C4", "#928EBE" , "lightgrey")

## Plot by site and date
ps.rare.top100@sam_data$Date_Long <- factor(ps.rare.top100@sam_data$Date_Long, levels = c("August 2018", "March 2019", "August 2019","March 2020"))

plot.relAbund_alv <- plot_bar(subset_samples(ps.rare.top100), x="CoralID", fill="Taxa") + 
  facet_grid(vars(Date_Long), vars(LTER_Site), scales = "free_x") + scale_fill_manual(values = colors)
plot.relAbund_alv <- plot.relAbund_alv + theme_bw(base_line_size = 1, base_rect_size = 1.5) + 
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 11.5), 
        axis.text.y = element_text(face = "bold", size = 11.5), title = element_text(face = "bold")) + 
  xlab("Coral ID") + ylab("Relative Abundance (%)") + labs(title = "Alveolates") + theme(strip.text = element_text(face = "bold", size = 11)) + labs(fill = "Taxa") 

## Plot by location and date
ps.rare.top100@sam_data$Location <- factor(ps.rare.top100@sam_data$Location, levels = c("North", "East", "West"))

plot.relAbund_alv_loc <- plot_bar(subset_samples(ps.rare.top100), x="CoralID", fill="Taxa") + 
  facet_grid(vars(Date_Long), vars(Location), scales = "free_x") + scale_fill_manual(values = colors)
plot.relAbund_alv_loc <- plot.relAbund_alv_loc + theme_bw(base_line_size = 1, base_rect_size = 1.5) + 
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 11.5), 
        axis.text.y = element_text(face = "bold", size = 11.5), title = element_text(face = "bold")) + 
  xlab("Coral ID") + ylab("Relative Abundance (%)") + labs(title = "Alveolates") + theme(strip.text = element_text(face = "bold", size = 11)) + labs(fill = "Taxa") 


plots <- ggarrange(plot.relAbund, plot.relAbund_alv, labels = c("A", "B"), ncol = 1, nrow = 2)
plots_loc <- ggarrange(plot.relAbund_loc, plot.relAbund_alv_loc, labels = c("A", "B"), ncol = 1, nrow = 2)

ggplot2::ggsave(here::here("plot_relabund.png"), plots,
                height = 600, width = 500, units = "mm",
                scale = 0.5, dpi = 1000)

ggplot2::ggsave(here::here("plot_relabund_loc.png"), plots_loc,
                height = 600, width = 500, units = "mm",
                scale = 0.5, dpi = 1000)


## Quantifying relative abundance of corallicola in different locations
### Use rarefied phyloseq object that was transformed to relative abundances

tax_glom(ps.rare.trans, taxrank = "Family", NArm = FALSE) -> ps.rare.trans.fam
ps.rare.cor <- subset_taxa(ps.rare.trans.fam, Family == "Corallicolidae")

melt <- psmelt(ps.rare.cor)
melt %>%
  select(Sample, Abundance, Location, Family, Date_Longer) -> melt

melt %>%
  filter(Date_Longer == "August 2018 (Before)") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = Location) %>%
  arrange()

melt %>%
  filter(Date_Longer == "March 2019 (During)") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = Location) %>%
  arrange()

melt %>%
  filter(Date_Longer == "August 2019 (After)") %>%
  summarise(mean = mean(Abundance), sd = sd(Abundance), .by = Location) %>%
  arrange()


kw_A18 <- broom::tidy(kruskal.test(Abundance ~ Location, data = subset(melt, Date_Longer == "August 2018 (Before)")))
kw_A18$Date <- "A18"

kw_M19 <- broom::tidy(kruskal.test(Abundance ~ Location, data = subset(melt, Date_Longer == "March 2019 (During)")))
kw_M19$Date <- "M19"
dunn_test(Abundance ~ Location, data = subset(melt, Date_Longer == "March 2019 (During)"), p.adjust.method = "fdr") # north and east sig diff

kw_A19 <- broom::tidy(kruskal.test(Abundance ~ Location, data = subset(melt, Date_Longer == "August 2019 (After)")))
kw_A19$Date <- "A19"
dunn_test(Abundance ~ Location, data = subset(melt, Date_Longer == "August 2019 (After)"), p.adjust.method = "fdr") # north and west sig diff 

kw_M20 <- broom::tidy(kruskal.test(Abundance ~ Location, data = subset(melt, Date_Longer == "March 2020 (1 Year Later)")))
kw_M20$Date <- "M20"

rbind(kw_A18, kw_M19, kw_A19, kw_M20) -> kw_output
write.csv(kw_output, here::here("kw_output.csv"))
