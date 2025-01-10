# Author: Sunni Patton
# Last edited: 11/8/24
# Title: Creating metadata file
# Overview: Built a metadata file based on original file name structure to use in phyloseq object

## Set seed ====
set.seed(123)

## Load libraries ====
library(ggplot2)
library(dplyr)
library(rstatix)
library(broom)

## Extract sample name data ====
### Read in seq_table_nometa 
readRDS(here::here("Output/00 - Read Preprocessing Output/seq_table_nometa.rds")) -> seq_table_nometa
sampleID <- rownames(seq_table_nometa)

sampleID_2 <- sampleID
sampleID_2 <- gsub("5_5", "5.5", sampleID_2)
sampleID_2 <- gsub("3_5", "3.5", sampleID_2)

site <- as.character(sapply(strsplit(sampleID_2, "_"), `[`,1))

reef <- as.character(sapply(strsplit(sampleID_2, "_"), `[`,2))

coralID <- as.character(sapply(strsplit(sampleID, "_ACR_"), `[`,2))
coralID <- as.character(sapply(strsplit(coralID, "_"), `[`,1))

date <- as.character(sapply(strsplit(sampleID, "_ACR_"), `[`,2))
date <- as.character(sapply(strsplit(date, "_"), `[`,2))

## Create metadata dataframe ====
metadata <- data.frame(SampleID = sampleID, LTER_Site = site, Reef_Type = reef, CoralID = coralID, Date = date)

metadata$Month[metadata$Date == "M20" |
                 metadata$Date == "M19"] <-"March"
metadata$Month[metadata$Date == "A18" |
                 metadata$Date == "A19"] <-"August"

metadata$Year[metadata$Date == "A18"] <- "2018"
metadata$Year[metadata$Date == "M19" | 
                metadata$Date == "A19"] <- "2019"
metadata$Year[metadata$Date == "M20"] <- "2020"

metadata$Date_Long <- paste0(metadata$Month, " ", metadata$Year)

metadata$Date_Longer[metadata$Date == "A18"] <- "August 2018 (Before)"
metadata$Date_Longer[metadata$Date == "M19"] <- "March 2019 (During)"
metadata$Date_Longer[metadata$Date == "A19"] <- "August 2019 (After)"
metadata$Date_Longer[metadata$Date == "M20"] <- "March 2020 (1 Year Later)"

metadata$Location[metadata$LTER_Site == "0" |
                    metadata$LTER_Site == "1" |
                    metadata$LTER_Site == "2"] <- "North"
metadata$Location[metadata$LTER_Site == "3" |
                    metadata$LTER_Site == "3.5" |
                    metadata$LTER_Site == "4"] <- "East"
metadata$Location[metadata$LTER_Site == "5" |
                    metadata$LTER_Site == "5.5"] <- "West"

write.csv(metadata, here::here("metadata.csv"))

