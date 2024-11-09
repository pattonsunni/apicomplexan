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

## Load data ====
readRDS(here::here("Output/02 - Phyloseq Preprocessing Output/ps.rare.rds")) -> ps.rare

## Agglomerate to genus level ====
ps.rare.genus <- tax_glom(ps.rare, taxrank = "Genus", NArm=FALSE) 

## Calculate alpha diversity ====
estimate_richness(ps.rare.genus, measures = c("Shannon", "Chao1")) -> alphaDiv

### Add Samples column
alphaDiv$SampleID <- rownames(ps.rare.genus@sam_data)

## Save alpha diversity output 
write.csv(alphaDiv, here::here("Output/05 - Alpha Diversity Output/alphaDiv.csv"))

### Combine alpha diversity with sam_data 
merge(as.matrix(ps.rare.genus@sam_data), alphaDiv, by = "SampleID") -> alphaDiv 

## Assess normality and variance assumptions ====
### Histogram
hist(alphaDiv$Shannon) # Not normally dist
### Normality
shapiro.test(alphaDiv$Shannon) # Not normally dist 
gofTest(alphaDiv$Shannon, distribution = "gamma") #appears to have gamma dist  

## Plot shannon and chao1 by site ====
alphaDiv$Date_Longer <- factor(alphaDiv$Date_Longer, c("August 2018 (Before)", "March 2019 (During)", "August 2019 (After)", "March 2020 (1 Year Later)"), ordered = TRUE)

plot_alphadiv_shan <-  ggplot(alphaDiv, aes(x = LTER_Site, y = Shannon, color = Date_Longer)) + 
  scale_color_manual(values = c("#1289ff", "#ff2c2c", "#ff8d0e", "#098943")) + geom_boxplot(lwd = 1.1, outlier.colour = "NA")
plot_alphadiv_shan <- plot_alphadiv_shan + geom_point(aes(color = Date_Longer), alpha = 0.5, 
                                            position = position_jitterdodge(jitter.width = 0.1)) + 
  theme_bw() + theme(axis.text = element_text(face = "bold", size = 11.5), 
                     axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold")) +
  labs(color = "Timepoint") + xlab("LTER Site")
  
plot_alphadiv_chao1 <-  ggplot(alphaDiv, aes(x = LTER_Site, y = log(Chao1), color = Date_Longer)) + 
  scale_color_manual(values = c("#1289ff", "#ff2c2c", "#ff8d0e", "#098943")) + geom_boxplot(lwd = 1.1, outlier.colour = "NA")
plot_alphadiv_chao1 <- plot_alphadiv_chao1 + geom_point(aes(color = Date_Longer), alpha = 0.5, 
                                                      position = position_jitterdodge(jitter.width = 0.1)) + 
  theme_bw() + theme(axis.text = element_text(face = "bold", size = 11.5), 
                     axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold")) +
  labs(color = "Timepoint") + xlab("LTER Site") + ylab("Log Transformed Chao1 Diversity")

## Linear models - Shannon ====
mod1 <- lmer(Chao1 ~  LTER_Site + (1|CoralID), data = alphaDiv)
anova(mod1) # p = 0.08 Shannon not significant by site

mod2 <- lmer(Shannon ~ Date + LTER_Site + Date*LTER_Site + (1|CoralID), data = alphaDiv)
shapiro.test(residuals(mod2)) # residuals normally distributed
anova(mod2) # Date and Date*LTER_Site significant 

### pairwise comparisons
mod_shan_pw <- emmeans(mod2, pairwise ~ LTER_Site | Date, adjust = "tukey") 
mod_shan_pw <- broom::tidy(mod_shan_pw$contrasts)
mod_shan_pw <- rstatix::add_significance(data = as.data.frame(mod_shan_pw), p.col = "adj.p.value")
mod_shan_pw$Metric <- "Shannon"

## Linear models - Chao1 ====
mod3 <- lmer(log(Chao1) ~  LTER_Site + (1|CoralID), data = alphaDiv)
anova(mod3) # p = 0.508 Chao1 not significant by site

mod4 <- lmer(log(Chao1) ~ Date + LTER_Site + Date*LTER_Site + (1|CoralID), data = alphaDiv)
shapiro.test(residuals(mod4)) # residuals normally distributed
anova(mod4) # Date*LTER_Site significant

### pairwise comparisons
mod_chao1_pw <- emmeans(mod2, pairwise ~ LTER_Site | Date, adjust = "tukey") 
mod_chao1_pw <- broom::tidy(mod_chao1_pw$contrasts)
mod_chao1_pw <- rstatix::add_significance(data = as.data.frame(mod_chao1_pw), p.col = "adj.p.value")
mod_chao1_pw$Metric <- "Chao1"

alpha_pairwise <- rbind(mod_shan_pw, mod_chao1_pw)
write.csv(alpha_pairwise, here::here("Output/04 - Alpha Diversity Output/alpha_pairwise.csv"))