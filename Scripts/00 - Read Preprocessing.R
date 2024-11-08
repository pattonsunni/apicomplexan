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
filtFs <- file.path(path, "filterAndTrim", paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filterAndTrim", paste0(sampleNames, "_R_filt.fastq.gz"))


## Quality filtering and trimming ====
### Based on quality profile, trim to c(270, 220)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(270, 220),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)
saveRDS(out, here::here("Output/00 - Read Preprocessing Output/out.rds"))

## Assess reads lost ====
sum(out[,1])-sum(out[,2]) #4,305,326  reads lost 

filter_fun = data.frame(out)
filter_fun
ratio = sum(filter_fun$reads.out)/sum(filter_fun$reads.in)
ratio # ~30% reads lost

## Learn errors ====
errF <- learnErrors(filtFs, multithread = FALSE) # 112757940  total bases in 417622  reads from 2 samples will be used for learning the error rates 
saveRDS(errF, here::here("Output/00 - Read Preprocessing Output/errF.rds"))
errR <- learnErrors(filtRs, multithread = FALSE) # 135082640 total bases in 614012 reads from 3 samples will be used for learning the error rates
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
dim(seq_table)# 35215 contigs in 89 samples

table(nchar(getSequences(seq_table))) # large range of contig length, but a majority fall at the correct length (396 - 407)

### Only keep those contigs that are from 399-407bp
seq_table <- seq_table[,nchar(colnames(seq_table)) %in% 396:407]
table(nchar(getSequences(seq_table)))
dim(seq_table) # In 89 samples, we have 32029 contigs
sum(seq_table) # 14,076,706 reads in these 32029 contigs

saveRDS(seq_table, here::here("Output/00 - Read Preprocessing Output/seq_table.rds"))

## Remove chimeric sequences from sequence table ====
seq_table_nochim <- removeBimeraDenovo(seq_table, method="consensus", multithread=TRUE, verbose=TRUE) # Identified 30760  bimeras out of 32029  input sequences
dim(seq_table_nochim) # in 89 samples, we now have 1269 contigs
sum(seq_table) - sum(seq_table_nochim) # 3,543,197 reads corresponding to 30760  contigs removed

## We still have 10,533,509 reads left

saveRDS(seq_table_nochim, here::here("Output/00 - Read Preprocessing Output/seq_table_nochim.rds"))

## Assign taxonomy and create taxonomy table ====
taxa <- assignTaxonomy(seq_table_nochim, here::here("pr2_version_5.0.0_SSU_dada2.fasta.gz"), multithread = TRUE, verbose = TRUE)
# multithread should be set to TRUE (even on windows) except when running filterAndTrim

saveRDS(taxa, here::here("Output/00 - Read Preprocessing Output/taxa.rds"))

## Change taxa column names to reflect correct taxonomy ====
colnames(taxa) <- c("Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species")

## Remove metazoan sequences ====
dim(seq_table_nochim) # 1,269 contigs
sum(seq_table_nochim) # 10533509 reads

# Make new sequence table 
is.meta <- taxa[,"Subdivision"] %in% "Metazoa"
seq_table_nometa <- seq_table_nochim[,!is.meta]
dim(seq_table_nometa) # Now we have 1222 sequence variants in our 89 samples (47 removed)
sum(seq_table_nochim) - sum(seq_table_nometa) # 47 contigs corresponding to 206,541 reads

# Make new taxonomy table
taxa_nometa <- taxa[!is.meta,]
dim(taxa_nometa)

## Remove bacterial sequences ====
# Make new sequence table
is.bact <- taxa_nometa[,"Domain"] %in% "Bacteria"
seq_table_nobact <- seq_table_nometa[,!is.bact]
dim(seq_table_nobact) # Now we have 1220 sequence variants in our 89 samples (2 removed)
sum(seq_table_nometa) - sum(seq_table_nobact) # 2 contigs corresponding to 49 reads
saveRDS(seq_table_nobact, here::here("Output/00 - Read Preprocessing Output/seq_table_nobact.rds"))

# Make new taxonomy table
taxa_nobact <- taxa_nometa[!is.bact,]
dim(taxa_nobact)
saveRDS(taxa_nobact, here::here("Output/00 - Read Preprocessing Output/taxa_nobact.rds"))

## Assess total number of reads removed and ASV frequency ====
sum(seq_table_nobact) # 10326919 reads left after quality control
dim(seq_table_nobact) # 1220 contigs left after quality control

summary(colSums(seq_table_nobact)) 
# Minimum number of times an ASV shows up in a sample is 1, 1st quartile: 8, median: 66.5, mean: 8464.7, 3rd quartile: 530, max: 2391129
summary(rowSums(seq_table_nobact)) # Minimum reads in a sample is 65, 1st quart: 90468 median: 117190 mean: 116033 3rd quant: 136847 max: 592625 

# See which sample had that low read number
sort(rowSums(seq_table_nobact)) # Sample with 65 reads is a negative control 

# Tracking reads removed at each step ====
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaForward, getN), sapply(dadaReverse, getN), sapply(contigs, getN), rowSums(seq_table_nochim), rowSums(seq_table_nometa), rowSums(seq_table_nobact))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim", "nometa", "nobact")

# Update sample names based on track rows (simply forcing column names to be sampleNames incorrectly assigns the names)
sampleNames_new <- sapply(strsplit(basename(rownames(track)), "-"), `[`,7) # Removes all beginning information, but still left with information at the end (S#_R1_001.fastq.gz)
sampleNames_new <- sapply(strsplit(basename(sampleNames_new), "\\."), `[`,1) # Removes fastq.gz at the end, but still left with the _S*_R1_001
sampleNames_new <- gsub("_S\\d+\\d?\\d?\\d?_R1_001", "", sampleNames_new) 
sampleNames_new

# Assign new row names (shorter sample name) using sampleNames_new
rownames(track) <- sampleNames_new
track$SampleID <- sampleNames_new

write.csv(track, here::here("Output/00 - Read Preprocessing Output/track.csv"))