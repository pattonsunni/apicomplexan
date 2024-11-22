# Author: Sunni Patton
# Last edited: 11/9/24
# Title: Beta Diversity
# Overview: Calculating beta diversity metrics

## Set seed ====
set.seed(123)

## Load libraries ====
library(microViz)
library(pairwiseAdonis)
library(phyloseq)
library(ggpubr)
library(cowplot)

## Load data ====
readRDS(here::here("Output/02 - Phyloseq Preprocessing Output/ps.All.rds")) -> ps.All

## Validate phyloseq object and fix if necessary ====
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

## Make ps object without forereef samples ====
ps.sub <- subset_samples(ps.All, Reef_Type != "FOR")

## Does LTER Site or Location affect beta diversity prior to disturbance (A18)? ====
## Subset data 
ps.A18 <- subset_samples(ps.sub, Date == "A18") # 24 total samples 
View(ps.A18@sam_data) 
# Site 0: 3 samples, site 1: 3 samples, site 2: 4 samples, site 3: 4 samples, site 3.5: 4 samples, site 4: 3 samples, site 5: 1 sample, site 5.5: 2 samples
# North: 10 samples, east: 11 samples, west: 3 samples

## rCLR transform data
ps.A18.rclr <- ps.A18 %>%
  tax_transform(trans = "rclr", rank = "Genus") %>%
  ord_calc(method = "PCA")

## Calculate Euclidean distance
A18_rclr_euc <- phyloseq::distance(ps.A18.rclr, method = "euclidean")

## Extract sample data
sampledf_A18 <- data.frame(sample_data(ps.A18.rclr))

## adonis2 (permanova)
adonis2_A18_site <- broom::tidy(adonis2(A18_rclr_euc ~ LTER_Site, data = sampledf_A18, perm = 999)) # p = 0.061; beta diversity is not significantly different by LTER Site
adonis2_A18_site$Comparison <- "LTER Site"
adonis2_A18_site$Date <- "August 2018 (Before)"

adonis2_A18_location <- broom::tidy(adonis2(A18_rclr_euc ~ Location, data = sampledf_A18, perm = 999)) # p = 0.026; beta diversity significantly different by Location
adonis2_A18_location$Comparison <- "Location"
adonis2_A18_location$Date <- "August 2018 (Before)"

## Pairwise adonis 
pwadonis_A18_loc <- pairwise.adonis(A18_rclr_euc, sampledf_A18$Location, p.adjust.m = "fdr", perm = 999) # North vs east significantly different
pwadonis_A18_loc$Date <- "August 2018 (Before)"


## Does LTER Site or Location affect beta diversity during disturbance (M19)? ====
## Subset data 
ps.M19 <- subset_samples(ps.sub, Date == "M19") # 25 total samples 
View(ps.M19@sam_data) 
# Site 0: 3 samples, site 1: 3 samples, site 2: 6 samples, site 3: 2 samples, site 3.5: 4 samples, site 4: 3 samples, site 5: 2 samples, site 5.5: 2 samples
# North: 12 samples, east: 9 samples, west: 4 samples

## rCLR transform data
ps.M19.rclr <- ps.M19 %>%
  tax_transform(trans = "rclr", rank = "Genus") %>%
  ord_calc(method = "PCA")

## Calculate Euclidean distance
M19_rclr_euc <- phyloseq::distance(ps.M19.rclr, method = "euclidean")

## Extract sample data
sampledf_M19 <- data.frame(sample_data(ps.M19.rclr))

## adonis2 (permanova)
adonis2_M19_site <- broom::tidy(adonis2(M19_rclr_euc ~ LTER_Site, data = sampledf_M19, perm = 999)) # p = 0.059; beta diversity is not significantly different by LTER Site
adonis2_M19_site$Comparison <- "LTER Site"
adonis2_M19_site$Date <- "March 2019 (During)"

adonis2_M19_location <- broom::tidy(adonis2(M19_rclr_euc ~ Location, data = sampledf_M19, perm = 999)) # p = 0.006; beta diversity significantly different by Location
adonis2_M19_location$Comparison <- "Location"
adonis2_M19_location$Date <- "March 2019 (During)"

## Pairwise adonis 
pwadonis_M19_loc <- pairwise.adonis(M19_rclr_euc, sampledf_M19$Location, p.adjust.m = "fdr", perm = 999) # North vs east significantly different
pwadonis_M19_loc$Date <- "March 2019 (During)"


## Does LTER Site or Location affect beta diversity after disturbance (A19)? ====
## Subset data 
ps.A19 <- subset_samples(ps.sub, Date == "A19") # 22 total samples 
View(ps.A19@sam_data) 
# Site 0: 1 samples, site 1: 3 samples, site 2: 5 samples, site 3: 2 samples, site 3.5: 4 samples, site 4: 3 samples, site 5: 2 samples, site 5.5: 2 samples
# North: 9 samples, east: 9 samples, west: 4 samples

## rCLR transform data
ps.A19.rclr <- ps.A19 %>%
  tax_transform(trans = "rclr", rank = "Genus") %>%
  ord_calc(method = "PCA")

## Calculate Euclidean distance
A19_rclr_euc <- phyloseq::distance(ps.A19.rclr, method = "euclidean")

## Extract sample data
sampledf_A19 <- data.frame(sample_data(ps.A19.rclr))

## adonis2 (permanova)
adonis2_A19_site <- broom::tidy(adonis2(A19_rclr_euc ~ LTER_Site, data = sampledf_A19, perm = 999)) # p = 0.032; beta diversity is significantly different by LTER Site
adonis2_A19_site$Comparison <- "LTER Site"
adonis2_A19_site$Date <- "August 2019 (After)"

adonis2_A19_location <- broom::tidy(adonis2(A19_rclr_euc ~ Location, data = sampledf_A19, perm = 999)) # p = 0.005; beta diversity significantly different by Location
adonis2_A19_location$Comparison <- "Location"
adonis2_A19_location$Date <- "August 2019 (After)"


## Pairwise adonis 
pwadonis_A19_site <- pairwise.adonis(A19_rclr_euc, sampledf_A19$LTER_Site, p.adjust.m = "fdr", perm = 999) # Not significant after p adj
pwadonis_A19_site$Date <- "August 2019 (After)"

pwadonis_A19_loc <- pairwise.adonis(A19_rclr_euc, sampledf_A19$Location, p.adjust.m = "fdr", perm = 999) # North vs east and north vs west significantly different
pwadonis_A19_loc$Date <- "August 2019 (After)"

## Does LTER Site or Location affect beta diversity one year after disturbance (M20)? ====
## Subset data 
ps.M20 <- subset_samples(ps.sub, Date == "M20") # 12 total samples 
View(ps.M20@sam_data) 
# Site 0: 0 samples, site 1: 1 samples, site 2: 4 samples, site 3: 0 samples, site 3.5: 2 samples, site 4: 2 samples, site 5: 2 samples, site 5.5: 1 samples
# North: 5 samples, east: 4 samples, west: 3 samples

## rCLR transform data
ps.M20.rclr <- ps.M20 %>%
  tax_transform(trans = "rclr", rank = "Genus") %>%
  ord_calc(method = "PCA")

## Calculate Euclidean distance
M20_rclr_euc <- phyloseq::distance(ps.M20.rclr, method = "euclidean")

## Extract sample data
sampledf_M20 <- data.frame(sample_data(ps.M20.rclr))

## adonis2 (permanova)
adonis2_M20_site <- broom::tidy(adonis2(M20_rclr_euc ~ LTER_Site, data = sampledf_M20, perm = 999)) # p = 0.177; beta diversity is not significantly different by LTER Site
adonis2_M20_site$Comparison <- "LTER Site"
adonis2_M20_site$Date <- "March 2020 (1 Year Later)"

adonis2_M20_location <- broom::tidy(adonis2(M20_rclr_euc ~ Location, data = sampledf_M20, perm = 999)) # p = 0.02; beta diversity significantly different by Location
adonis2_M20_location$Comparison <- "Location"
adonis2_M20_location$Date <- "March 2020 (1 Year Later)"

## Pairwise adonis 
pwadonis_M20_loc <- pairwise.adonis(M20_rclr_euc, sampledf_M20$Location, p.adjust.m = "fdr", perm = 999) # Nothing significant after p adjustment 
pwadonis_M20_loc$Date <- "March 2020 (1 Year Later)"

## Combine statistics output ====
adonis2_output_date <- rbind(adonis2_A18_site, adonis2_A18_location, adonis2_M19_site, adonis2_M19_location, adonis2_A19_site, adonis2_A19_location, adonis2_M20_site, adonis2_M20_location)
write.csv(adonis2_output_date, here::here("Output/05 - Beta Diversity Output/adonis_output_date_final_sub.csv"))

pwadonis_output_date <- rbind(pwadonis_A18_loc, pwadonis_M19_loc, pwadonis_A19_loc, pwadonis_A19_site, pwadonis_M20_loc)
write.csv(pwadonis_output_date, here::here("Output/05 - Beta Diversity Output/pwadonis_output_date_final_sub.csv"))

## Determine the affect of time on beta diversity in the north shore ====
## Subset data
ps.north <- subset_samples(ps.All, Location == "North") # 36 samples total
View(ps.north@sam_data) # A18: 10 samples, M19: 12 samples, A19: 9 samples, M20: 5 samples

## rCLR transform data
ps.north.rclr <- ps.north %>%
  tax_transform(trans = "rclr", rank = "Genus") %>%
  ord_calc(method = "PCA")

## Calculate Euclidean distance
north_rclr_euc <- phyloseq::distance(ps.north.rclr, method = "euclidean")

## Extract sample data
sampledf_north <- data.frame(sample_data(ps.north.rclr))

## adonis2 (permanova)
adonis2_north_site <- broom::tidy(adonis2(north_rclr_euc ~ LTER_Site, data = sampledf_north, perm = 999)) # p = 0.076; within the north shore, no significant differences observed in beta diversity by site
adonis2_north_site$Comparison <- "LTER Site"
adonis2_north_site$Location <- "North"

adonis2_north_date <- broom::tidy(adonis2(north_rclr_euc ~ Date, data = sampledf_north, perm = 999)) # p = 0.001; beta diversity significantly different by date
adonis2_north_date$Comparison <- "Date"
adonis2_north_date$Location <- "North"

## Pairwise adonis 
pw_north_date <- pairwise.adonis(north_rclr_euc, sampledf_north$Date, p.adjust.m = "fdr", perm = 999) # sig between during vs after (p = 0.015), before vs after (p = 0.006), and before vs 1 year after (p = 0.026)
pw_north_date$Location <- "North"

## Determine the affect of time on beta diversity in the east shore ====
## Subset data
ps.east <- subset_samples(ps.All, Location == "East") # 33 samples total
View(ps.east@sam_data) # A18: 11 samples, M19: 9 samples, A19: 9 samples, M20: 4 samples

## rCLR transform data
ps.east.rclr <- ps.east %>%
  tax_transform(trans = "rclr", rank = "Genus") %>%
  ord_calc(method = "PCA")

## Calculate Euclidean distance
east_rclr_euc <- phyloseq::distance(ps.east.rclr, method = "euclidean")

## Extract sample data
sampledf_east <- data.frame(sample_data(ps.east.rclr))

## adonis2 (permanova)
adonis2_east_site <- broom::tidy(adonis2(east_rclr_euc ~ LTER_Site, data = sampledf_east, perm = 999)) # p = 0.091; within the east shore, no significant differences observed in beta diversity by site
adonis2_east_site$Comparison <- "LTER Site"
adonis2_east_site$Location <- "East"

adonis2_east_date <- broom::tidy(adonis2(east_rclr_euc ~ Date, data = sampledf_east, perm = 999)) # p = 0.384; beta diversity not significantly different by date in east shore
adonis2_east_date$Comparison <- "Date"
adonis2_east_date$Location <- "East"

## Determine the affect of time on beta diversity in the west shore ====
## Subset data
ps.west <- subset_samples(ps.All, Location == "West") # 14 samples total
View(ps.west@sam_data) # A18: 3 samples, M19: 4 samples, A19: 4 samples, M20: 3 samples

## rCLR transform data
ps.west.rclr <- ps.west %>%
  tax_transform(trans = "rclr", rank = "Genus") %>%
  ord_calc(method = "PCA")

## Calculate Euclidean distance
west_rclr_euc <- phyloseq::distance(ps.west.rclr, method = "euclidean")

## Extract sample data
sampledf_west <- data.frame(sample_data(ps.west.rclr))

## adonis2 (permanova)
adonis2_west_site <- broom::tidy(adonis2(west_rclr_euc ~ LTER_Site, data = sampledf_west, perm = 999)) # p = 0.49; within the west shore, no significant differences observed in beta diversity by site
adonis2_west_site$Comparison <- "LTER Site"
adonis2_west_site$Location <- "West"

adonis2_west_date <- broom::tidy(adonis2(west_rclr_euc ~ Date, data = sampledf_west, perm = 999)) # p = 0.452; beta diversity not significantly different by date in east shore
adonis2_west_date$Comparison <- "Date"
adonis2_west_date$Location <- "West"

## Determine the affect of time on beta diversity in the east and west shores combined ====
## Subset data
ps.eastwest <- subset_samples(ps.All, Location != "North") # 47 samples total
View(ps.eastwest@sam_data) # A18: 14 samples, M19: 13 samples, A19: 13 samples, M20: 7 samples

## rCLR transform data
ps.eastwest.rclr <- ps.eastwest %>%
  tax_transform(trans = "rclr", rank = "Genus") %>%
  ord_calc(method = "PCA")

## Calculate Euclidean distance
eastwest_rclr_euc <- phyloseq::distance(ps.eastwest.rclr, method = "euclidean")

## Extract sample data
sampledf_eastwest <- data.frame(sample_data(ps.eastwest.rclr))

## adonis2 (permanova)
adonis2_eastwest_site <- broom::tidy(adonis2(eastwest_rclr_euc ~ LTER_Site, data = sampledf_eastwest, perm = 999)) # p = 0.253; within the east and west shores combined, no significant differences observed in beta diversity by site
adonis2_eastwest_site$Comparison <- "LTER Site"
adonis2_eastwest_site$Location <- "East + West"

adonis2_eastwest_date <- broom::tidy(adonis2(eastwest_rclr_euc ~ Date, data = sampledf_eastwest, perm = 999)) # p = 0.184; beta diversity not significantly different by date in east shore
adonis2_eastwest_date$Comparison <- "Date"
adonis2_eastwest_date$Location <- "East + West"


## Combine statistics output ====
## adonis2 
adonis2_output <- rbind(adonis2_north_site, adonis2_north_date, adonis2_east_site, adonis2_east_date, adonis2_west_site, adonis2_west_date)
write.csv(adonis2_output, here::here("Output/05 - Beta Diversity Output/adonis_output_loc_sub_final.csv"))

write.csv(pw_north_date, here::here("Output/05 - Beta Diversity Output/pwadonis_output_loc_sub_final.csv"))

## Plot ====

## North shore over time
ps.north.rclr@sam_data$Date_Longer <- factor(ps.north.rclr@sam_data$Date_Longer, c("August 2018 (Before)", "March 2019 (During)", "August 2019 (After)", "March 2020 (1 Year Later)", 
                                                 order = TRUE))
PCA_north <- ps.north.rclr %>%
  ord_plot( 
    colour = "Date_Longer",
    plot_taxa = FALSE, 
    auto_caption = NA
  ) + stat_ellipse(aes(group = Date_Longer, color = Date_Longer), linewidth = 1.5) + theme_bw() + scale_color_brewer(palette = "Dark2") +
  theme_bw(base_line_size = 1.5, base_rect_size = 1) + 
  theme(axis.text = element_text(face = "bold", size = 14), axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) + theme(strip.text = element_text(face = "bold", size = 14)) + guides(color = guide_legend(title = "Date")) + theme(legend.position = "none") + facet_wrap(~Location)
  # + theme(legend.position = c(0.8, 0.08)) + theme(legend.key.size = unit(1, "mm")) +

## East shore over time
ps.east.rclr@sam_data$Date_Longer <- factor(ps.east.rclr@sam_data$Date_Longer, c("August 2018 (Before)", "March 2019 (During)", "August 2019 (After)", "March 2020 (1 Year Later)", 
                                                                                           order = TRUE))
PCA_east <- ps.east.rclr %>%
  ord_plot( 
    colour = "Date_Longer",
    plot_taxa = FALSE, 
    auto_caption = NA
  ) + stat_ellipse(aes(group = Date_Longer, color = Date_Longer), linewidth = 1.5) + theme_bw() + scale_color_brewer(palette = "Dark2") +
  theme_bw(base_line_size = 1.5, base_rect_size = 1) + 
  theme(axis.text = element_text(face = "bold", size = 14), axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) + facet_wrap(~Location) + theme(strip.text = element_text(face = "bold", size = 14)) + guides(color = guide_legend(title = "Date")) + theme(legend.position = "none") 
  #  + theme(legend.position = c(0.8, 0.06)) + theme(legend.key.size = unit(1, "mm"))

## West shore over time
ps.west.rclr@sam_data$Date_Longer <- factor(ps.west.rclr@sam_data$Date_Longer, c("August 2018 (Before)", "March 2019 (During)", "August 2019 (After)", "March 2020 (1 Year Later)", 
                                                                                         order = TRUE))
PCA_west <- ps.west.rclr %>%
  ord_plot( 
    colour = "Date_Longer",
    plot_taxa = FALSE, 
    auto_caption = NA
  ) + stat_ellipse(aes(group = Date_Longer, color = Date_Longer), linewidth = 1.5) + theme_bw() + scale_color_brewer(palette = "Dark2") +
  theme_bw(base_line_size = 1.5, base_rect_size = 1) + 
  theme(axis.text = element_text(face = "bold", size = 14), axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) + facet_wrap(~Location) + theme(strip.text = element_text(face = "bold", size = 14)) + guides(color = guide_legend(title = "Date")) + theme(legend.position = "none") 
  # + theme(legend.position = c(0.8, 0.06)) + theme(legend.key.size = unit(1, "mm"))

## Supplemental plot
ggarrange(PCA_north, PCA_east, PCA_west, ncol = 3) -> betadiv_shores_plot

ggplot2::ggsave(here::here("Output/05 - Beta Diversity Output/01 - beta_div.png"), betadiv_plot,
                height = 450, width = 800, units = "mm",
                scale = 0.5, dpi = 1000)

## Plot north, and east+west
ps.north.rclr@sam_data$Date_Longer <- factor(ps.north.rclr@sam_data$Date_Longer, c("August 2018 (Before)", "March 2019 (During)", "August 2019 (After)", "March 2020 (1 Year Later)", 
                                                                                   order = TRUE))
PCA_north <- ps.north.rclr %>%
  ord_plot( 
    colour = "Date_Longer",
    plot_taxa = FALSE, 
    auto_caption = NA
  ) + stat_ellipse(aes(group = Date_Longer, color = Date_Longer), linewidth = 1.5) + theme_bw() + scale_color_brewer(palette = "Dark2") +
  theme_bw(base_line_size = 1.5, base_rect_size = 1) + 
  theme(axis.text = element_text(face = "bold", size = 14), axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) + theme(strip.text = element_text(face = "bold", size = 14)) +
  guides(color = guide_legend(title = "Date")) + theme(legend.position = c(0.8, 0.08)) + theme(legend.key.size = unit(1, "mm")) + labs(title = "North Shore")


ps.eastwest.rclr@sam_data$Date_Longer <- factor(ps.eastwest.rclr@sam_data$Date_Longer, c("August 2018 (Before)", "March 2019 (During)", "August 2019 (After)", "March 2020 (1 Year Later)", 
                                                                                         order = TRUE))
PCA_eastwest <- ps.eastwest.rclr %>%
  ord_plot( 
    colour = "Date_Longer",
    plot_taxa = FALSE, 
    auto_caption = NA
  ) + stat_ellipse(aes(group = Date_Longer, color = Date_Longer), linewidth = 1.5) + theme_bw() + scale_color_brewer(palette = "Dark2") +
  theme_bw(base_line_size = 1.5, base_rect_size = 1) + 
  theme(axis.text = element_text(face = "bold", size = 14), axis.title = element_text(face = "bold", size = 14), 
        title = element_text(face = "bold")) + theme(strip.text = element_text(face = "bold", size = 14, title)) +
  guides(color = guide_legend(title = "Date")) + theme(legend.position = c(0.8, 0.08)) + theme(legend.key.size = unit(1, "mm")) + labs(title = "East and West Shores")

ggarrange(PCA_north, PCA_eastwest) -> betadiv_plot2
ggplot2::ggsave(here::here("Output/05 - Beta Diversity Output/02 - beta_div2.png"), betadiv_plot2,
                height = 400, width = 600, units = "mm",
                scale = 0.5, dpi = 1000)














## Beta dispersion ====
## North shore
disp_north <- betadisper(north_rclr_euc, sampledf_north$Date, bias.adjust = TRUE, type = "centroid")
df_disp_north <- as.data.frame(disp_north$distances)
df_disp_north$Samples <- paste0(rownames(df_disp_north))

date <- as.character(sapply(strsplit(df_disp_north$Samples, "_ACR_"), `[`,2))
date <- as.character(sapply(strsplit(date, "_"), `[`,2))

df_disp_north$Date <- date

df_disp_north$Date[df_disp_north$Date == "A18"] <-"August 2018 (Before)"
df_disp_north$Date[df_disp_north$Date == "M19"] <-"March 2019 (During)"
df_disp_north$Date[df_disp_north$Date == "A19"] <-"August 2019 (After)"
df_disp_north$Date[df_disp_north$Date == "M20"] <-"March 2020 (1 Year Later)"

colnames(df_disp_north) <- c("Distances", "Samples", "Date")
df_disp_north$Location <- "North"


broom::tidy(anova(disp_north)) # Dispersion not significant across time in north shore (p = 0.542)
north_stats <- broom::tidy(p.adjust(permutest(betadisper(north_rclr_euc, sampledf_north$Date, type = "centroid", 
                              bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
north_stats <- add_significance(north_stats, p.col = "x")

north_stats$Location <- "North"
colnames(north_stats) <- c("Comparison", "p.adj", "p.adj.signif", "Location")
north_stats$group1 <- sapply(strsplit(north_stats$Comparison, "-"), `[`,1)
north_stats$group2 <- sapply(strsplit(north_stats$Comparison, "-"), `[`,2)

north_stats$y.position <- 16

north_stats$group1[north_stats$group1 == "A18"] <-"August 2018 (Before)"
north_stats$group1[north_stats$group1 == "M19"] <-"March 2019 (During)"
north_stats$group1[north_stats$group1 == "A19"] <-"August 2019 (After)"
north_stats$group1[north_stats$group1 == "M20"] <-"March 2020 (1 Year Later)"

north_stats$group2[north_stats$group2 == "A18"] <-"August 2018 (Before)"
north_stats$group2[north_stats$group2 == "M19"] <-"March 2019 (During)"
north_stats$group2[north_stats$group2 == "A19"] <-"August 2019 (After)"
north_stats$group2[north_stats$group2 == "M20"] <-"March 2020 (1 Year Later)"


## East shore
### LTER Site
disp_east <- betadisper(east_rclr_euc, sampledf_east$Date, bias.adjust = TRUE, type = "centroid")
df_disp_east <- as.data.frame(disp_east$distances)
df_disp_east$Samples <- paste0(rownames(df_disp_east))

date <- as.character(sapply(strsplit(df_disp_east$Samples, "_ACR_"), `[`,2))
date <- as.character(sapply(strsplit(date, "_"), `[`,2))

df_disp_east$Date <- date

df_disp_east$Date[df_disp_east$Date == "A18"] <-"August 2018 (Before)"
df_disp_east$Date[df_disp_east$Date == "M19"] <-"March 2019 (During)"
df_disp_east$Date[df_disp_east$Date == "A19"] <-"August 2019 (After)"
df_disp_east$Date[df_disp_east$Date == "M20"] <-"March 2020 (1 Year Later)"

colnames(df_disp_east) <- c("Distances", "Samples", "Date")
df_disp_east$Location <- "East"


broom::tidy(anova(disp_east)) # Dispersion significant over time in east shore (p = 0.002714)
east_stats <- broom::tidy(p.adjust(permutest(betadisper(east_rclr_euc, sampledf_east$Date, type = "centroid", 
                              bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr')) # Dispersion significant between A18-A19 (p = 0.003), A19-M19 (p = 0.003)

east_stats <- add_significance(east_stats, p.col = "x")

east_stats$Location <- "East"
colnames(east_stats) <- c("Comparison", "p.adj", "p.adj.signif", "Location")
east_stats$group1 <- sapply(strsplit(east_stats$Comparison, "-"), `[`,1)
east_stats$group2 <- sapply(strsplit(east_stats$Comparison, "-"), `[`,2)

east_stats$y.position <- 16

east_stats$group1[east_stats$group1 == "A18"] <-"August 2018 (Before)"
east_stats$group1[east_stats$group1 == "M19"] <-"March 2019 (During)"
east_stats$group1[east_stats$group1 == "A19"] <-"August 2019 (After)"
east_stats$group1[east_stats$group1 == "M20"] <-"March 2020 (1 Year Later)"

east_stats$group2[east_stats$group2 == "A18"] <-"August 2018 (Before)"
east_stats$group2[east_stats$group2 == "M19"] <-"March 2019 (During)"
east_stats$group2[east_stats$group2 == "A19"] <-"August 2019 (After)"
east_stats$group2[east_stats$group2 == "M20"] <-"March 2020 (1 Year Later)"


## West shore
disp_west <- betadisper(west_rclr_euc, sampledf_west$Date, bias.adjust = TRUE, type = "centroid")
df_disp_west <- as.data.frame(disp_west$distances)
df_disp_west$Samples <- paste0(rownames(df_disp_west))

date <- as.character(sapply(strsplit(df_disp_west$Samples, "_ACR_"), `[`,2))
date <- as.character(sapply(strsplit(date, "_"), `[`,2))

df_disp_west$Date <- date

df_disp_west$Date[df_disp_west$Date == "A18"] <-"August 2018 (Before)"
df_disp_west$Date[df_disp_west$Date == "M19"] <-"March 2019 (During)"
df_disp_west$Date[df_disp_west$Date == "A19"] <-"August 2019 (After)"
df_disp_west$Date[df_disp_west$Date == "M20"] <-"March 2020 (1 Year Later)"

colnames(df_disp_west) <- c("Distances", "Samples", "Date")
df_disp_west$Location <- "West"

anova(disp_west) # beta diversity not significant over time in west shore
west_stats <- broom::tidy(p.adjust(permutest(betadisper(west_rclr_euc, sampledf_west$Date, type = "centroid", 
                                                        bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr')) # Dispersion significant between A18-A19 (p = 0.003), A19-M19 (p = 0.003)

west_stats <- add_significance(west_stats, p.col = "x")

west_stats$Location <- "West"
colnames(west_stats) <- c("Comparison", "p.adj", "p.adj.signif", "Location")
west_stats$group1 <- sapply(strsplit(west_stats$Comparison, "-"), `[`,1)
west_stats$group2 <- sapply(strsplit(west_stats$Comparison, "-"), `[`,2)

west_stats$y.position <- 16

west_stats$group1[west_stats$group1 == "A18"] <-"August 2018 (Before)"
west_stats$group1[west_stats$group1 == "M19"] <-"March 2019 (During)"
west_stats$group1[west_stats$group1 == "A19"] <-"August 2019 (After)"
west_stats$group1[west_stats$group1 == "M20"] <-"March 2020 (1 Year Later)"

west_stats$group2[west_stats$group2 == "A18"] <-"August 2018 (Before)"
west_stats$group2[west_stats$group2 == "M19"] <-"March 2019 (During)"
west_stats$group2[west_stats$group2 == "A19"] <-"August 2019 (After)"
west_stats$group2[west_stats$group2 == "M20"] <-"March 2020 (1 Year Later)"


## Combine dataframes ====
rbind(df_disp_north, df_disp_east, df_disp_west) -> disp_df

rbind(north_stats, east_stats, west_stats) -> disp_stats_df

disp_stats_df[7,]$y.position <- 16
disp_stats_df[10,]$y.position <- 18


## Plot beta dispersion ====
## set order
disp_df$Location <- factor(disp_df$Location, levels = c("North", "East", "West"), 
                            order = TRUE)

disp_df$Date <- factor(disp_df$Date, c("August 2018 (Before)", "March 2019 (During)", "August 2019 (After)", "March 2020 (1 Year Later)", 
                           order = TRUE))

plot_disp <-  ggplot(disp_df, aes(x = Date, y = Distances, color = Date))  + 
  geom_boxplot(lwd = 1.25, outlier.colour = "NA") + theme_bw(base_line_size = 1.5, base_rect_size = 1.75)

plot_disp <- plot_disp + geom_point(aes(color = Date), alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1)) +
   ylab("Distance to Centroid") + scale_color_brewer(palette = "Dark2")
plot_disp <- plot_disp + theme(axis.text = element_text(face = "bold", size = 14), 
                               axis.title = element_text(face = "bold", size = 14), 
                               title = element_text(face = "bold"), axis.title.x = element_blank(), 
                               axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  guides(color = guide_legend(title = "Date")) + theme(strip.text = element_text(face = "bold", size = 14)) +
  stat_pvalue_manual(disp_stats_df, label = "p.adj.signif", hide.ns = TRUE, size = 6) + 
  theme(legend.position = c(0.85, 0.9)) + theme(legend.key.size = unit(1, "mm")) + facet_grid(~factor(Location, levels=c('North', 'East', 'West'))) 
  


ggplot2::ggsave(here::here("Output/05 - Beta Diversity Output/02 - beta_disp.png"), plot_disp,
                height = 400, width = 600, units = "mm",
                scale = 0.5, dpi = 1000)

ggarrange(PCA_north, PCA_east, PCA_west, nrow = 1, ncol = 3) -> x

ggarrange(x, plot_disp, ncol = 1, nrow = 2, common.legend = TRUE, legend = "bottom") -> PCA_disp

ggplot2::ggsave(here::here("Output/05 - Beta Diversity Output/03 - PCA_disp.png"), PCA_disp,
                height = 400, width = 600, units = "mm",
                scale = 0.5, dpi = 1000)

## Save files ====
write.csv(disp_df, here::here("Output/05 - Beta Diversity Output/beta_disp_final.csv"))
write.csv(disp_stats_df, here::here("Output/05 - Beta Diversity Output/beta_disp_stats_final.csv"))
