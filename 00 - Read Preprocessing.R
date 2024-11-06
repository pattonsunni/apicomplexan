# Author: Sunni Patton
# Date: 02/12/24 (last edited)
# Title: Read preprocessing

# All reads have been run through cutadapt to remove forward, reverse, reverse complement forward, reverse complement reverse primers

## Load libraries ====
library(here)
library(dada2)
library(dplyr)

## Set seed for reproducibility ====
set.seed(123)

## Set path to files ====
path <- here::here("Data")

### Set path for forward and reverse reads
fnFs <- sort(list.files(path, pattern = "R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "R2_001.fastq.gz", full.names = TRUE))

## Extract sample names from files ====
### Assume name follows: lane1-s0001-index--Findex-Rindex-SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(basename(fnFs), "-"), `[`,7) # Removes all beginning information, but still left with information at the end (S#_R1_001.fastq.gz)
sampleNames <- sapply(strsplit(basename(sampleNames), "\\."), `[`,1) # Removes fastq.gz at the end, but still left with the _S*_R1_001
sampleNames <- gsub("_S\\d+\\d?\\d?\\d?_R1_001", "", sampleNames) 
sampleNames

## Set file destination for files after quality filtering ====
filtFs <- file.path(path, "filterAndTrim", paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filterAndTrim", paste0(sampleNames, "_R_filt.fastq.gz"))


## Quality filtering and trimming ====
### Based on quality profile, trim to c(270, 220)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(270, 220),
                        maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                        compress=TRUE, multithread=FALSE)
saveRDS(out, here::here("Output/00 - Read Preprocessing Output/out.rds"))

## if these reads below don't match, it was likely ran with 275,275 -> looks like I actually ran it like this

## Assess reads lost ====
sum(out[,1])-sum(out[,2]) # 9,022,331 reads lost 

filter_fun = data.frame(out)
filter_fun
ratio = sum(filter_fun$reads.out)/sum(filter_fun$reads.in)
ratio # ~50% reads lost 

## Learn errors ====
errF <- learnErrors(filtFs, multithread = FALSE) # 124737250 total bases in 453590 reads from 3 samples will be used for learning the error rates 
saveRDS(errF, here::here("Output/00 - Read Preprocessing Output/errF.rds"))
errR <- learnErrors(filtRs, multithread = FALSE) # 124737250 total bases in 453590 reads from 3 samples will be used for learning the error rates
saveRDS(errR, here::here("Output/00 - Read Preprocessing Output/errR.rds"))

## Ensure sample naming is consistent ====
names(filtFs)<-sampleNames
names(filtRs)<-sampleNames

## Infer sample sequence (dada) ====
dadaForward <- dada(filtFs, err=errF, multithread=FALSE)
saveRDS(dadaForward, here::here("Output/00 - Read Preprocessing Output/dadaForward.rds"))
dadaReverse <- dada(filtRs, err=errR, multithread=FALSE)
saveRDS(dadaReverse, here::here("Output/00 - Read Preprocessing Output/dadaReverse.rds"))

## Create contigs and sequence table ====
contigs <- mergePairs(dadaForward, filtFs, dadaReverse, filtRs)
saveRDS(contigs, here::here("Output/00 - Read Preprocessing Output/contigs.rds"))

### Make sequence table and visualize contig length and frequency
seq_table <- makeSequenceTable(contigs) 
dim(seq_table) # Total of 4761 contigs


