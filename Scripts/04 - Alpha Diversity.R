# Author: Sunni Patton
# Last edited: 11/08/2024
# Title: Alpha Diversity
# Overview: Calculating alpha diversity

## Set seed ====
set.seed(123)

## Load libraries ====
library(phyloseq)
library(stats)
library(ggpubr)
library(rstatix)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(fitdistrplus)
library(EnvStats)
library(broom)
library(svglite)

## Load data ====
readRDS(here::here("Output/02 - Phyloseq Preprocessing Output/ps.rare.rds")) -> ps.rare

## Agglomerate to genus level ====
ps.rare.genus <- tax_glom(ps.rare, taxrank = "Genus", NArm=FALSE) 

## Calculate alpha diversity ====
estimate_richness(ps.rare.genus, measures = c("Shannon", "Chao1")) -> alphaDiv

### Add Samples column
alphaDiv$SampleID <- rownames(ps.rare.genus@sam_data)

### Combine alpha diversity with sam_data 
merge(as.matrix(ps.rare.genus@sam_data), alphaDiv, by = "SampleID") -> alphaDiv 

## Save alpha diversity output 
write.csv(alphaDiv, here::here("alphaDiv.csv"))


## Assess normality and variance assumptions ====
### Histogram
hist(alphaDiv$Shannon) # Not normally dist
### Normality
shapiro.test(alphaDiv$Shannon) # Not normally dist 
gofTest(alphaDiv$Shannon, distribution = "gamma") #appears to have gamma dist 

## Linear models - Shannon ====
mod1 <- lmer(Shannon ~ Date + LTER_Site + Date*LTER_Site + (1|CoralID), data = alphaDiv)
shapiro.test(residuals(mod1)) # residuals normally distributed
lmem_shan <- broom::tidy(anova(mod1)) # Date and Date*LTER_Site significant 
lmem_shan$Metric <- "Shannon"

### pairwise comparisons
mod_shan_pw <- emmeans(mod1, pairwise ~ LTER_Site | Date, adjust = "tukey") 
mod_shan_pw <- broom::tidy(mod_shan_pw$contrasts)
mod_shan_pw <- rstatix::add_significance(data = as.data.frame(mod_shan_pw), p.col = "adj.p.value")
mod_shan_pw$Metric <- "Shannon"

mod_shan_pw2 <- emmeans(mod1, pairwise ~ Date | LTER_Site, adjust = "tukey")
mod_shan_pw2 <- broom::tidy(mod_shan_pw2$contrasts)
mod_shan_pw2 <- rstatix::add_significance(data = as.data.frame(mod_shan_pw2), p.col = "adj.p.value")
mod_shan_pw2$Metric <- "Shannon"

## Linear models - Chao1 ====
mod2 <- lmer(log(Chao1) ~ Date + LTER_Site + Date*LTER_Site + (1|CoralID), data = alphaDiv)
shapiro.test(residuals(mod2)) # residuals normally distributed
lmem_chao1 <- broom::tidy(anova(mod2)) # Date*LTER_Site significant
lmem_chao1$Metric <- "Chao1"

### pairwise comparisons
mod_chao1_pw <- emmeans(mod2, pairwise ~ LTER_Site | Date, adjust = "tukey") 
mod_chao1_pw <- broom::tidy(mod_chao1_pw$contrasts)
mod_chao1_pw <- rstatix::add_significance(data = as.data.frame(mod_chao1_pw), p.col = "adj.p.value")
mod_chao1_pw$Metric <- "Chao1"

mod_chao1_pw2 <- emmeans(mod2, pairwise ~ Date | LTER_Site, adjust = "tukey")
mod_chao1_pw2 <- broom::tidy(mod_chao1_pw2$contrasts)
mod_chao1_pw2 <- rstatix::add_significance(data = as.data.frame(mod_chao1_pw2), p.col = "adj.p.value")
mod_chao1_pw2$Metric <- "Chao1"

lmem_all <- rbind(lmem_shan, lmem_chao1)
write.csv(lmem_all, here::here("lmem_output.csv"))

alpha_pairwise <- rbind(mod_shan_pw, mod_chao1_pw)
write.csv(alpha_pairwise, here::here("alpha_pairwise.csv"))

alpha_pairwise2 <- rbind(mod_shan_pw2, mod_chao1_pw2)
write.csv(alpha_pairwise2, here::here("alpha_pairwise2.csv"))









## Plot shannon and chao1 by site ====
### Prep statistics for plot
alpha_pairwise$group1 <- paste0(alpha_pairwise$contrast)
alpha_pairwise$group2 <- paste0(alpha_pairwise$contrast)

alpha_pairwise$group1 <- sapply(strsplit(basename(alpha_pairwise$group1), "-"), `[`,1)
alpha_pairwise$group1  <- gsub("LTER_Site", "", alpha_pairwise$group1) 

alpha_pairwise$group2 <- sapply(strsplit(basename(alpha_pairwise$group2), "-"), `[`,2)
alpha_pairwise$group2  <- gsub("LTER_Site", "", alpha_pairwise$group2) 

alpha_pairwise$y.position <- 2.1

colnames(alpha_pairwise)[10] <- "p.adj.signif"

as.numeric(alpha_pairwise$group1) -> alpha_pairwise$group1
as.numeric(alpha_pairwise$group2) -> alpha_pairwise$group2

as.factor(alphaDiv$LTER_Site) -> alphaDiv$LTER_Site
as.factor(alpha_pairwise$group1) -> alpha_pairwise$group1
as.factor(alpha_pairwise$group2) -> alpha_pairwise$group2

alpha_pairwise$y.position[31] <- 1
alpha_pairwise$y.position[37] <- 1.5
alpha_pairwise$y.position[42] <- 1.3
alpha_pairwise$y.position[47] <- 1.7
alpha_pairwise$y.position[48] <- 1.8
alpha_pairwise$y.position[49] <- 2
alpha_pairwise$y.position[50] <- 1.2

alpha_pairwise <- subset(alpha_pairwise, p.adj.signif != "ns" & adj.p.value != "NA")


alphaDiv$Date_Longer <- factor(alphaDiv$Date_Longer, c("August 2018 (Before)", "March 2019 (During)", "August 2019 (After)", "March 2020 (1 Year Later)"), ordered = TRUE)

plot_alphadiv_shan <-  ggplot(alphaDiv, aes(x = LTER_Site, y = Shannon, color = Date_Longer)) + 
  scale_color_brewer(palette = "Dark2") + geom_boxplot(lwd = 1.1, outlier.colour = "NA") 

plot_alphadiv_shan <- plot_alphadiv_shan + geom_point(aes(color = Date_Longer), size = 1.5, alpha = 0.5, 
                                            position = position_jitterdodge(jitter.width = 0.1)) + 
  theme_bw(base_line_size = 1.5, base_rect_size = 1.75) + theme(axis.text = element_text(face = "bold", size = 11.5), 
                     axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold")) +
  stat_pvalue_manual(subset(alpha_pairwise, Metric == "Shannon"), label = "p.adj.signif", hide.ns = TRUE, size = 6) +
  labs(color = "Timepoint") + xlab("LTER Site") 




alpha_pairwise$y.position[7] <- 4.5
alpha_pairwise$y.position[8] <- 4.3
alpha_pairwise$y.position[9] <- 3.7
alpha_pairwise$y.position[10] <- 3.8
alpha_pairwise$y.position[11] <- 3.5
alpha_pairwise$y.position[12] <- 3.3
alpha_pairwise$y.position[13] <- 3


plot_alphadiv_chao1 <-  ggplot(alphaDiv, aes(x = LTER_Site, y = log(Chao1), color = Date_Longer)) + 
  scale_color_brewer(palette = "Dark2") + geom_boxplot(lwd = 1.1, outlier.colour = "NA")
plot_alphadiv_chao1 <- plot_alphadiv_chao1 + geom_point(aes(color = Date_Longer), alpha = 0.5, 
                                                      position = position_jitterdodge(jitter.width = 0.1)) + 
  theme_bw(base_line_size = 1.5, base_rect_size = 1.75) + theme(axis.text = element_text(face = "bold", size = 11.5), 
                     axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold")) +
  stat_pvalue_manual(subset(alpha_pairwise, Metric == "Chao1"), label = "p.adj.signif", hide.ns = TRUE, size = 6) +
  labs(color = "Timepoint") + xlab("LTER Site") + ylab("Log Transformed Chao1 Diversity")


alpha_plot <- ggarrange(plot_alphadiv_shan, plot_alphadiv_chao1, ncol = 2, nrow = 1, labels = c("A", "B"), common.legend = TRUE)

ggplot2::ggsave(here::here("plot_alpha.svg"), alpha_plot,
                height = 300, width = 500, units = "mm",
                scale = 0.5, dpi = 1000)
 