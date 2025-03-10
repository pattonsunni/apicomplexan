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
library(dplyr)
library(effsize)

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


## Fix extra space around taxa names
as.data.frame(ps.rare@tax_table) -> ps.rare.taxtable
ps.rare.taxtable$Taxa[ps.rare.taxtable$Taxa == "Cladocopium  "] <- "Cladocopium"
ps.rare.taxtable$Taxa[ps.rare.taxtable$Taxa == "Cladocopium "] <- "Cladocopium"
ps.rare.taxtable$Taxa[ps.rare.taxtable$Taxa == "Symbiodinium "] <- "Symbiodinium"
ps.rare.taxtable$Taxa[ps.rare.taxtable$Taxa == "Corallicola aquarius "] <- "Corallicola aquarius"
ps.rare.taxtable$Taxa[ps.rare.taxtable$Taxa == "Rhodomelaceae X"] <- "Rhodomelaceae Family"

as.matrix(ps.rare.taxtable) -> ps.rare.taxtable

tax_table(ps.rare) <- ps.rare.taxtable


## Transform to relative abundance ====
ps.rare.trans <- transform_sample_counts(ps.rare, function(OTU) OTU/sum(OTU))

## Prune most abundant taxa ====
top100 <- names(sort(taxa_sums(ps.rare.trans), decreasing = TRUE))[1:100] 
ps.rare.top100 <- prune_taxa(top100, ps.rare.trans)

# Change coralID to factor
ps.rare.top100@sam_data$CoralID <- as.factor(ps.rare@sam_data$CoralID)

## Plot relative abundance for all non-metazoan micro-eukaryotes ====
tax_glom(ps.rare.top100, taxrank = "Taxa", NArm = FALSE) -> ps.rare.top100

colors <- c("#4BACC6", "#EA6312", "#800000",  "#ADB17D", "#FFB547", "#7F5F52", "#D092A7",
            "#9B6BF2", "#E33D6F", "#8EC0C1", "#4472C4", "#928EBE", "#24693D")



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
ps.rare.top100@sam_data$CoralID <- as.factor(ps.rare.top100@sam_data$CoralID)

## Plot relative abundance for alveolates ====
tax_glom(ps.rare.top100, taxrank = "Taxa", NArm = FALSE) -> ps.rare.top100

colors <- c("#800000","#A1D68B", "#ADB17D", "#FFB547", "#7F5F52", "#B31166", "#24693D", "#F3E79A", "#4472C4", "#928EBE" , "lightgrey")

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

ggplot2::ggsave(here::here("Output/03 - Relative Abundance Output/plot_relabund.png"), plots,
                height = 600, width = 650, units = "mm",
                scale = 0.5, dpi = 1000)

ggplot2::ggsave(here::here("Output/03 - Relative Abundance Output/plot_relabund_loc.png"), plots_loc,
                height = 600, width = 650, units = "mm",
                scale = 0.5, dpi = 1000)


## Quantifying relative abundance of corallicola in different locations ====
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
dunn_M19 <- as.data.frame(dunn_test(Abundance ~ Location, data = subset(melt, Date_Longer == "March 2019 (During)"), p.adjust.method = "fdr")) # north and east sig diff
dunn_M19$Date <- "M19"

kw_A19 <- broom::tidy(kruskal.test(Abundance ~ Location, data = subset(melt, Date_Longer == "August 2019 (After)")))
kw_A19$Date <- "A19"
dunn_A19 <- as.data.frame(dunn_test(Abundance ~ Location, data = subset(melt, Date_Longer == "August 2019 (After)"), p.adjust.method = "fdr")) # north and west sig diff 
dunn_A19$Date <- "A19"

kw_M20 <- broom::tidy(kruskal.test(Abundance ~ Location, data = subset(melt, Date_Longer == "March 2020 (1 Year Later)")))
kw_M20$Date <- "M20"

rbind(kw_A18, kw_M19, kw_A19, kw_M20) -> kw_output
write.csv(kw_output, here::here("Output/03 - Relative Abundance Output/kw_output.csv"))

rbind(dunn_M19, dunn_A19) -> dunn_output
write.csv(dunn_output, here::here("Output/03 - Relative Abundance Output/dunn_output.csv"))

## Effect size east vs north 
melt %>%
  filter(Date_Longer == "March 2019 (During)" & Location == "East")
melt %>%
  filter(Date_Longer == "March 2019 (During)" & Location == "North")

group1 <- c(0.23667723, 0.13854277, 0.13445970, 0.13206617, 0.08546287, 0.05040479, 0.03428370,0.02752552, 0.02731433)
group2 <- c(0.11587469, 0.10116156, 0.03583245, 0.03484688, 0.03097501, 0.02949666, 0.02604717, 0.02541359, 0.02316086, 0.02316086, 0.02231609 , 0.00499824)

cohen.d(group1, group2)

## Effect size west vs north 
melt %>%
  filter(Date_Longer == "March 2019 (During)" & Location == "West")

group1 <- c(0.12439282, 0.09123548, 0.08497008, 0.05251672)
group2 <- c(0.11587469, 0.10116156, 0.03583245, 0.03484688, 0.03097501, 0.02949666, 0.02604717, 0.02541359, 0.02316086, 0.02316086, 0.02231609 , 0.00499824)

cohen.d(group1, group2)


## Assessing corallicolidae total abundance ====
# check if corallicolidae are in all samples in both rarefied and unrarefied data
subset_taxa(ps.All, Family == "Corallicolidae") -> ps.All.cor
sort(sample_sums(ps.All.cor))
sum(taxa_sums(ps.All.cor)) # 644,291 reads are corallicolidae

subset_taxa(ps.rare, Family == "Corallicolidae") -> ps.rare.cor
sort(sample_sums(ps.rare.cor))
sum(taxa_sums(ps.rare.cor)) # 83,950 reads are corallicolidae

## Corallicolidae are found in all samples across the island in the nonrarfied data, although 4 samples have <100 reads
## After rarefying, corallicolidae are found in 82/83 samples, but now 14 samples have <100 reads


## Use nonrarefied data to plot total corallicolidae abundance by shore type over time
ps.All.cor.fam <- ps.All.cor %>%
  tax_glom("Family")

ps.All.cor.fam@sam_data$Location <- factor(ps.All.cor.fam@sam_data$Location, c("North", "East", "West"), ordered = TRUE)
ps.All.cor.fam@sam_data$Date_Long <- factor(ps.All.cor.fam@sam_data$Date_Long, c("August 2018", "March 2019", "August 2019", "March 2020"), ordered = TRUE)

cor_nonrare <- phyloseq::psmelt(ps.All.cor.fam) %>%
  ggplot(data = ., aes(x = Date_Long, y = Abundance, colour = Date_Long)) +
  geom_boxplot(lwd = 1.25, stat = "boxplot", outlier.color = "red", alpha = 0.5) +
  geom_jitter(height = 0, width = 0.2) + 
  facet_wrap(~Location) +
  theme_bw(base_line_size = 1.5, base_rect_size = 1.75) +
  scale_color_brewer(palette = "Dark2") + 
  theme(axis.text.x = element_text(angle = 50, face = "bold", size = 11.5, margin = margin(t = 20, l = 20)),
        axis.text.y = element_text(face = "bold", size = 11.5), title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 12), legend.position = "none") + 
  xlab("Date") + labs(title = "Nonrarefied") 

ps.All.cor.fam@sam_data$Date_Long <- factor(ps.All.cor.fam@sam_data$Date_Long, c("August 2018", "March 2019", "August 2019", "March 2020"), ordered = TRUE)

phyloseq::psmelt(ps.All.cor.fam) %>%
  ggplot(data = ., aes(x = Location, y = Abundance, colour = Location)) +
  geom_boxplot(lwd = 1.25, stat = "boxplot", outlier.color = "red", alpha = 0.5) +
  geom_jitter(height = 0, width = 0.2) + 
  facet_wrap(~Date_Long, ncol = 4, nrow = 1) +
  theme_bw(base_line_size = 1.5, base_rect_size = 1.75) +
  scale_color_brewer(palette = "Dark2") + 
  theme(axis.text.x = element_text(face = "bold", size = 11.5),
        axis.text.y = element_text(face = "bold", size = 11.5), title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 12), legend.position = "none") + 
  xlab("Date") + labs(title = "Nonrarefied")

## Use rarefied data
ps.rare.cor.fam <- ps.rare.cor %>%
  tax_glom("Family")

ps.rare.cor.fam@sam_data$Location <- factor(ps.rare.cor.fam@sam_data$Location, c("North", "East", "West"), ordered = TRUE)
ps.rare.cor.fam@sam_data$Date_Long <- factor(ps.rare.cor.fam@sam_data$Date_Long, c("August 2018", "March 2019", "August 2019", "March 2020"), ordered = TRUE)

cor_rare <- phyloseq::psmelt(ps.rare.cor.fam) %>%
  ggplot(data = ., aes(x = Date_Long, y = Abundance, colour = Date_Long)) +
  geom_boxplot(lwd = 1.25, stat = "boxplot", outlier.color = "red", alpha = 0.5) +
  geom_jitter(height = 0, width = 0.2) + 
  facet_wrap(~Location) +
  theme_bw(base_line_size = 1.5, base_rect_size = 1.75) +
  scale_color_brewer(palette = "Dark2") + 
  theme(axis.text.x = element_text(angle = 50, face = "bold", size = 11.5, margin = margin(t = 20, l = 20)),
        axis.text.y = element_text(face = "bold", size = 11.5), title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 12), legend.position = "none") + 
  xlab("Date") + labs(title = "Rarefied") 

ggarrange(cor_nonrare, cor_rare) -> cor_abund

ggplot2::ggsave(here::here("Output/03 - Relative Abundance Output/plot_cora_abund.svg"), cor_abund,
                height = 400, width = 700, units = "mm",
                scale = 0.5, dpi = 1000)
