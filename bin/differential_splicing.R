## internal information requires R/3.2.3 on cluster
## Rscript --vanilla [file_name] [working_dir] [project_name] [meta_data] [num_cores]

library(DEXSeq)
library(BiocParallel)

# get command line arguments
args <- commandArgs(TRUE)
print(args)

# set working directory
working_dir <- args[1]
setwd(working_dir)

# get project name
project_name <- args[2]

meta_data <-args[3]

num_cores <- as.numeric(args[4])

# load data
countData <- read.csv(paste(working_dir, "/", project_name, "_sj_counts", sep=""), header=TRUE, row.names=1, sep="\t",check.names=F)
print("countData")
dim(countData)
head(countData)
# >head(countData)
#                         112c_JM 488e_JM 501b_JM 635b_JM CEMM_05_CC1 CEMM_12_TK
# 20:49551774-49552684:-      52     126     149     136         108         99
# 20:49552800-49557401:-      51     100     153     125         176        110
# 20:49552800-49558567:-       4       2       0       2           8          6
# 20:49557446-49558567:-       0       0       0       0           0          0
# 20:49557471-49557641:-       2       2       7       4          15          0
# 20:49557471-49557665:-       9      12       8      23          28          2

featureIDs <- read.csv(
  paste(working_dir, "/", project_name, "_featureIDs", sep=""), 
  header=TRUE, sep="\t"
)
print("featureIDs")
head(featureIDs)
# > head(featureIDs)
# exon_id
# 1       1
# 2       2
# 3       3
# 4       4
# 5       5
# 6       6
groupIDs <- read.csv(
  paste(working_dir, "/", project_name, "_groupID", sep=""), 
  header=TRUE, sep="\t"
)
print("groupIDs")
head(groupIDs)
# > head(groupIDs)
# gene_id
# 1 ENSG00000000419
# 2 ENSG00000000419
# 3 ENSG00000000419
# 4 ENSG00000000419
# 5 ENSG00000000419
# 6 ENSG00000000419

transcripts <- read.csv(
  paste(working_dir, "/", project_name, "_transcripts", sep=""), 
  header=TRUE, sep="\t"
)
print("transcripts")
head(transcripts)
# > head(transcripts)
# X            transcripts
# 1 0 20:49551774-49552684:-
#   2 1 20:49552800-49557401:-
#   3 2 20:49552800-49558567:-
#   4 3 20:49557446-49558567:-
#   5 4 20:49557471-49557641:-
#   6 5 20:49557471-49557665:-
#   
sampleData <- read.csv(meta_data, header=TRUE, row.names=1, sep="\t")
print("sampleData")
dim(sampleData)
head(sampleData)

# some data transformation
featureIDs$exon_ID <- paste("exon", featureIDs$exon_id, sep="")
featureIDs <- as.vector(featureIDs$exon_ID)
groupIDs <- as.vector(groupIDs[,])
# #featureRanges <- read.csv(paste(dir, "exonIntervals", sep=""), header=TRUE)
transcripts <- as.vector(transcripts[,])

# reorder sample data so that it has the same order as columns in countData
sampleData.sorted <- sampleData[colnames(countData),]
colnames(sampleData.sorted) <- c("condition", "batch_name", "disease_abbr", "genotype")

# print reformatted data
print("Reformatted data:")
head(featureIDs)
head(groupIDs)
head(transcripts)
head(sampleData.sorted)

## DEXSeq workflow

print("Creating DEXSeq data object:") 
ECS = DEXSeqDataSet(
  countData=countData,
  sampleData=sampleData.sorted,
  design=formula(~ sample + exon + condition:exon),
  featureID=featureIDs,
  groupID=groupIDs,
  transcripts=transcripts)

# Pre-filtering
# dds <- dds[ rowSums(counts(dds)) > 1, ]

print("Estimating size factors:")
ECS = estimateSizeFactors(ECS)

print("Estimate Dispersion:")
BPPARAM=MulticoreParam(workers=num_cores)
ECS = estimateDispersions(ECS, BPPARAM=BPPARAM)
#ECS = fitDispersionFunction(ECS)
print("Test for differential expression:")
ECS = testForDEU(ECS, BPPARAM=BPPARAM)
print("Estimate exon fold changes:")
ECS = estimateExonFoldChanges(ECS, fitExpToVar = "condition")

print("Save R session:")
save.image(paste(paste(working_dir, "/", project_name, ".RData", sep="")))

print("Export results to file:")
counts <- counts(ECS, normalized=FALSE)
write.csv(counts, paste(working_dir, "/", project_name, "_counts", sep=""), sep="\t")

normalized_counts <- counts(ECS, normalized=TRUE)
write.csv(normalized_counts, paste(working_dir, "/", project_name, "_counts_normalized", sep=""), sep="\t")

dexseq_results <- DEXSeqResults(ECS)
write.csv(dexseq_results, paste(working_dir, "/", project_name, "_dexseq_results", sep=""), sep="\t")

## TODO
# vst function from DESeq2
# MDS

