# Author: Sunni Patton
# Date: 11/08/24 (last edited)
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
filtFs <- file.path(path, "filterAndTrim_121324", paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filterAndTrim_121324", paste0(sampleNames, "_R_filt.fastq.gz"))


## Quality filtering and trimming ====
### Based on quality profile, trim to c(270, 220)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(270, 220),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)
saveRDS(out, here::here("out.rds"))

## Assess reads lost ====
sum(out[,1])-sum(out[,2]) #19774760 - 15469434 -> 4305326 reads lost

filter_fun = data.frame(out)
filter_fun
ratio = sum(filter_fun$reads.out)/sum(filter_fun$reads.in)
ratio # ~32% reads lost

## Learn errors ====
errF <- learnErrors(filtFs, multithread = FALSE) # 112757940   total bases in 417622   reads from 2 samples will be used for learning the error rates 
saveRDS(errF, here::here("errF.rds"))
errR <- learnErrors(filtRs, multithread = FALSE) # 135082640  total bases in 614012  reads from 3 samples will be used for learning the error rates
saveRDS(errR, here::here("errR.rds"))


## Ensure sample naming is consistent ====
names(filtFs)<-sampleNames
names(filtRs)<-sampleNames

## Infer sample sequence (dada) ====
dadaForward <- dada(filtFs, err=errF, multithread=FALSE)
saveRDS(dadaForward, here::here("dadaForward.rds"))
dadaReverse <- dada(filtRs, err=errR, multithread=FALSE)
saveRDS(dadaReverse, here::here("dadaReverse.rds"))

## Create contigs and sequence table ====
contigs <- mergePairs(dadaForward, filtFs, dadaReverse, filtRs)
saveRDS(contigs, here::here("contigs.rds"))

### Make sequence table and visualize contig length and frequency
seq_table <- makeSequenceTable(contigs) 
dim(seq_table)# 35215 contigs in 89 samples

table(nchar(getSequences(seq_table))) # large range of contig length, but a majority fall at the correct length (396 - 406)

### Only keep those contigs that are from 399-406bp
seq_table <- seq_table[,nchar(colnames(seq_table)) %in% 396:406]
table(nchar(getSequences(seq_table)))
dim(seq_table) # In 89 samples, we have 31968 contigs
sum(seq_table) # 14075736 reads in these 31968 contigs

saveRDS(seq_table, here::here("seq_table.rds"))

## Remove chimeric sequences from sequence table ====
seq_table_nochim <- removeBimeraDenovo(seq_table, method="consensus", multithread=TRUE, verbose=TRUE) # Identified 30705    bimeras out of 31968  input sequences
dim(seq_table_nochim) # in 89 samples, we now have 1263 contigs
sum(seq_table) - sum(seq_table_nochim) # 3542736 reads corresponding to 31968 contigs removed

## We still have 10533000 reads left

saveRDS(seq_table_nochim, here::here("seq_table_nochim.rds"))

## Assign taxonomy and create taxonomy table ====
taxa <- assignTaxonomy(seq_table_nochim, here::here("pr2_version_5.0.0_SSU_dada2.fasta.gz"), multithread = TRUE, verbose = TRUE)
# multithread should be set to TRUE (even on windows) except when running filterAndTrim

saveRDS(taxa, here::here("taxa.rds"))

dim(taxa) # 1263 ASVs
sum(seq_table_nochim) # 10533000 reads

## Remove bacterial sequences ====
is.bact <- taxa[,"Kingdom"] %in% "Bacteria"
seq_table_nobact <- seq_table_nochim[,!is.bact]

dim(seq_table_nobact) # 1261 ASVs
sum(seq_table_nobact) # 10532951 reads (49 reads were bacteria)

taxa_nobact <- taxa[!is.bact,]
dim(taxa_nobact) #1261 taxa

## Change taxa column names to reflect correct taxonomy ====
colnames(taxa_nobact) <- c("Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species")

## Remove metazoan sequences ====
# Make new sequence table 
is.meta <- taxa_nobact[,"Subdivision"] %in% "Metazoa"
seq_table_nometa <- seq_table_nobact[,!is.meta]
dim(seq_table_nometa) # Now we have 1216 sequence variants in our 89 samples (45 ASVs removed)
sum(seq_table_nobact) - sum(seq_table_nometa) # 45 contigs corresponding to 206484 reads

# Make new taxonomy table
taxa_nometa <- taxa_nobact[!is.meta,]
dim(taxa_nometa)

dim(seq_table_nometa) # 1216 contigs
sum(seq_table_nometa) # 10326467 reads

saveRDS(taxa_nometa, here::here("taxa_nometa.rds"))
saveRDS(seq_table_nometa, here::here("seq_table_nometa.rds"))

## Assess total number of reads removed and ASV frequency ====
summary(colSums(seq_table_nometa)) 
# Minimum number of times an ASV shows up in a sample is 1, 1st quartile: 8, median: 66.5, mean: 8492.2, 3rd quartile: 536.8 , max: 2391129.0  
summary(rowSums(seq_table_nometa)) # Minimum reads in a sample is 65, 1st quart: 90468     median: 117188     mean: 116028   3rd quant: 136845     max: 592613  

# See which sample had that low read number
sort(rowSums(seq_table_nometa)) # Sample with 65 reads is a negative control 

# Tracking reads removed at each step ====
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaForward, getN), sapply(dadaReverse, getN), sapply(contigs, getN), rowSums(seq_table_nochim), rowSums(seq_table_nobact),rowSums(seq_table_nometa))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim", "nobact", "nometa")

# Update sample names based on track rows (simply forcing column names to be sampleNames incorrectly assigns the names)
sampleNames_new <- sapply(strsplit(basename(rownames(track)), "-"), `[`,7) # Removes all beginning information, but still left with information at the end (S#_R1_001.fastq.gz)
sampleNames_new <- sapply(strsplit(basename(sampleNames_new), "\\."), `[`,1) # Removes fastq.gz at the end, but still left with the _S*_R1_001
sampleNames_new <- gsub("_S\\d+\\d?\\d?\\d?_R1_001", "", sampleNames_new) 
sampleNames_new

# Assign new row names (shorter sample name) using sampleNames_new
rownames(track) <- sampleNames_new
as.data.frame(track) -> track
track$SampleID <- sampleNames_new

write.csv(track, here::here("track.csv"))
