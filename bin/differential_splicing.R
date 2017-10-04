#! /usr/bin/env Rscript
#
#SBATCH --job-name=psi
#SBATCH --cpus-per-task=32
#SBATCH --distribution=block
#SBATCH --mem=256000
#SBATCH --time=14-00:00:00
#SBATCH --partition=longq
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --account=fschischlik
#SBATCH -o psi.%j.out
#SBATCH -e psi.%j.err
#
# R script to run differential splicing analysis
#
# requires R/3.2.3 on cluster
# Rscript --vanilla 
# [meta_data]
# [exp_name]
# [design_formula]
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
reformat_sampleData <- FALSE
verbose <- TRUE

# set working directory
setwd(working_dir)

ref.sampleData <- function(sampleData){
  
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
  
  return(sampleData)
}

if(!exists("ECS")){
  
  ## countData
  ############
  countData <- read.csv(
    paste(working_dir, "/", 
          project_name, "_sj_counts", 
          sep=""), 
    header=TRUE, 
    row.names=1, 
    sep="\t",
    check.names=F)
  
  ## featureIDs
  #############
  featureIDs <- read.csv(
    paste(working_dir, "/", 
          project_name, "_featureIDs", 
          sep=""), 
    header=TRUE, 
    sep="\t")
  
  featureIDs$exon_ID <- paste(
    "exon", 
    featureIDs$exon_id, 
    sep="")
  featureIDs <- as.vector(featureIDs$exon_ID)
  
  ## groupIDs
  ###########
  groupIDs <- read.csv(
    paste(working_dir, "/", 
          project_name, "_groupID", 
          sep=""), 
    header=TRUE, 
    sep="\t")
  
  groupIDs <- as.vector(groupIDs[,])
  
  ## transcripts
  ##############
  transcripts <- read.csv(
    paste(working_dir, "/", 
          project_name, "_transcripts", 
          sep=""), 
    header=TRUE, 
    sep="\t")
  
  transcripts <- as.vector(transcripts[,])
  
  ## sampleData
  #############
  sampleData <- read.csv(
    meta_data, 
    header=TRUE, 
    row.names=1, 
    sep="\t")
  
  # reorder sample data so that it has the same 
  # order as columns in countData
  sampleData.sorted <- sampleData[colnames(countData),]
  
  if (reformat_sampleData){
    sampleData <- ref.sampleData(sampleData)
  }
  
  ## print loaded resource files
  ##############################
  print("countData: ")
  dim(countData)
  head(countData)
  
  print("featureIDs: ")
  head(featureIDs)
  
  print("groupIDs")
  head(groupIDs)
  
  print("transcripts")
  head(transcripts)
  
  print("sampleData")
  dim(sampleData.sorted)
  head(sampleData.sorted)
  
  ## read in design formulas
  ##########################
  design <- read.csv(
    design_formulas,
    header=TRUE, 
    sep="\t")
  
  fitExpToVar <- as.character(
    design[design$exp_name==exp_name,]$fitExpToVar)
  formulaFullModel <- as.formula(
    as.character(
      design[design$exp_name==exp_name,]$formulaFullModel))
  formulaReducedModel <- as.formula(
    as.character(
      design[design$exp_name==exp_name,]$formulaReducedModel))
  
  print("Formulas:")
  print(formulaFullModel)
  print(formulaReducedModel)
  
  print("FitExpToVar:")
  print(fitExpToVar)
  
  if(verbose){
    
    ## itermediate save
    ###################
    print("Save R session before DEXSeq data object creation:")
    save.image(paste(
      working_dir, "/", 
      project_name, "_", 
      exp_name,"_","loaded_resources", 
      ".RData", 
      sep=""))
  }
  
  ## DEXSeq workflow
  ##################
  BPPARAM=MulticoreParam(workers=num_cores)
  
  print("Creating DEXSeq data object:") 
  ECS = DEXSeqDataSet(
    countData=countData,
    sampleData=sampleData.sorted,
    design=formulaFullModel,
    featureID=featureIDs,
    groupID=groupIDs,
    transcripts=transcripts)
  
  if(exp_name=="interaction1"){
    print(formulaFullModel)
    print(colData(ECS))
    m1 <- model.matrix(formulaFullModel, colData(ECS))
    all.zero <- apply(m1, 2, function(x) all(x==0))
    idx <- which(all.zero)
    m1 <- m1[,-idx]
    formulaFullModel <- m1
    ECS = DEXSeqDataSet(
      countData=countData,
      sampleData=sampleData.sorted,
      design=formulaFullModel,
      featureID=featureIDs,
      groupID=groupIDs,
      transcripts=transcripts)
  }
  
  ## free up some memory
  ######################
  rm(countData,
     featureIDs,
     groupIDs,
     transcripts)
  
  # Pre-filtering
  ###############
  # dds <- dds[ rowSums(counts(dds)) > 1, ]
  
  print("Estimating size factors:")
  ECS = estimateSizeFactors(ECS)
  
  print("Estimate Dispersion:")
  ECS = estimateDispersions(
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
  ECS = estimateExonFoldChanges(
    ECS, 
    fitExpToVar = fitExpToVar)
  
  ## save results
  ###############
  print("Save R session:")
  save.image(paste(paste(working_dir, "/", 
                         project_name, "_", 
                         exp_name, ".RData", sep="")))
  
  print("Export results to file:")
  counts <-featureCounts(
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
  
  normalized_counts <-featureCounts(
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
  
}

## continue with plotting if R object available
###############################################
if(!exists("ECS")){
  load(paste(paste(working_dir, "/", 
                   project_name, "_", 
                   exp_name, ".RData", sep="")))
}