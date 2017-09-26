## internal information requires R/3.2.3 on cluster
## Rscript --vanilla 
# [this_script] 
# [meta_data]
# [exp_name]
# [design_formulas]
# [working_dir] 
# [project_name] 
# [num_cores]

library(DEXSeq)
library(BiocParallel)
library(matrixStats)

# get command line arguments
args <- commandArgs(TRUE)
print(args)

meta_data <-args[1]
exp_name <- args[2]
design_formulas <- args[3]
working_dir <- args[4]
project_name <- args[5]
num_cores <- as.numeric(args[6])

# set working directory
setwd(working_dir)

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

sampleData <- read.csv(
  meta_data, 
  header=TRUE, 
  row.names=1, sep="\t")

print("sampleData")
dim(sampleData)
head(sampleData)

sampleData$bleeding_date <- as.Date(
  as.character(sampleData$bleeding_date),
  format="%d.%m.%Y")
sampleData$extraction_date <- as.Date(
  as.character(sampleData$extraction_date),
  format="%d.%m.%Y")
sampleData$date_of_birth <- as.Date(
  as.character(sampleData$date_of_birth),
  format="%d.%m.%Y")
sampleData$storage_time <- sampleData$bleeding_date - sampleData$extraction_date
sampleData$storage_time <- difftime(
  sampleData$extraction_date, 
  sampleData$bleeding_date, unit="days")
sampleData$storage_time <- as.numeric(
  sampleData$storage_time)
sampleData$age <- difftime(
  sampleData$bleeding_date, 
  sampleData$date_of_birth, unit="weeks")
sampleData$age <- as.numeric(sampleData$age)/52.14

#some data transformation
featureIDs$exon_ID <- paste("exon", featureIDs$exon_id, sep="")
featureIDs <- as.vector(featureIDs$exon_ID)
groupIDs <- as.vector(groupIDs[,])
# #featureRanges <- read.csv(paste(dir, "exonIntervals", sep=""), header=TRUE)
transcripts <- as.vector(transcripts[,])

# reorder sample data so that it has the same order as columns in countData
sampleData.sorted <- sampleData[colnames(countData),]

# read in design formulas
design <- read.csv(design_formulas)

fitExpToVar <- design[design$exp_name==exp_name,]$fitExpToVar
formulaFullModel <- design[design$exp_name==exp_name,]$formulaFullModel
formulaReducedModel <- design[design$exp_name==exp_name,]$formulaReducedModel

# print reformatted data
print("Reformatted data:")
head(featureIDs)
head(groupIDs)
head(transcripts)
head(sampleData.sorted)

## DEXSeq workflow

# assign number of cores
BPPARAM=MulticoreParam(workers=num_cores)

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
dxd_e1=estimateDispersions(
  ECS,
  formula=formulaFullModel,
  BPPARAM=BPPARAM)

print("Test for differential expression:")
ECS = testForDEU(
  ECS, 
  reducedModel=formulaReducedModel,
  fullModel=formulaFullModel,
  BPPARAM=BPPARAM)

print("Estimate exon fold changes:")
ECS = estimateExonFoldChanges(ECS, fitExpToVar = fitExpToVar)

print("Save R session:")
save.image(paste(paste(working_dir, "/", 
                       project_name, "_", 
                       exp_name, ".RData", sep="")))

print("Export results to file:")
dexseq_counts <-featureCounts(
  ECS, 
  normalized=FALSE)
write.table(
  counts, 
  file=paste(working_dir, "/", 
             project_name, "_",
             exp_name, "_counts", sep=""),
  quote=FALSE,
  col.names=TRUE,
  sep="\t")

dexseq_counts <-featureCounts(
  ECS, 
  normalized=TRUE)
write.table(
  normalized_counts, 
  file=paste(working_dir, "/", 
             project_name, "_",
             exp_name, "_counts_normalized", sep=""),
  quote=FALSE,
  col.names=TRUE,
  sep="\t")

dexseq_results <- DEXSeqResults(ECS)
write.table(
  dexseq_results,
  paste(working_dir, "/", 
        project_name, "_",
        exp_name, "_dexseq_results", sep=""),
  quote=FALSE,
  col.names=TRUE,
  sep="\t")
