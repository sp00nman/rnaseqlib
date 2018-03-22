#' Annotate variants from vcf2table format
#'
#'@param variants A dataframe with variants to filter \code{variants}
#'@param filter.samples A vector with name of healthy samples\code{filter.samples}
#'@param alt.count Filter criteria, minimum number of alternate reads \code{alt.count}
#'@param qualitybydepth Filter criteria, minimum quality by depth for variants \code{qualitybydepth}
#'@param vaf.healthy Minimum variant allele frequency required for healthy samples \code{vaf.healthy}
#'@param verbose For trouble shooting \code {verbose}
#'
#'@return vector with annotations, comma separated
#'
#'@keywords keywords
#'
#'@export
#'
#'@examples
#'
#'annotateVariants(rcm, filter.samples=c("C1", "C2"), alt.count=3, vaf=0.4, verbose=FALSE)
#'

annotateVariants <- function(
  variants, 
  filter.samples=healthy.samples,
  alt.count=2,
  qualitybydepth=2,
  vaf.healthy=0.4,
  verbose=FALSE){
  
  # these filters are the same as used for the patient mutation matrix  
  
  .Annotate <-function(row_data){
    if(verbose){print(row_data[["X1000g2015feb_eur"]])}
    if(verbose){print(str(row_data[["X1000g2015feb_eur"]]))}
    
    filter.annotation <- c()
    if(row_data[["PLOTNAME"]] %in% as.character(healthy.samples)){
      filter.annotation <- cbind(filter.annotation, "HEALTHY")
    }
    if(row_data[["CANONICAL_TRANSCRIPT_AAChange"]] %in% healthy.variants){
      filter.annotation <- cbind(filter.annotation, "INHEALTHY") 
    }
    if(row_data[["X1000g2015feb_eur"]] > 0.01){
      filter.annotation <- cbind(filter.annotation, "1000gEUR")
    }
    if(row_data[["X1000g2015feb_amr"]] > 0.01){
      filter.annotation <- cbind(filter.annotation, "1000gAMR")
    }
    if(row_data[["X1000g2015feb_eas"]] > 0.01){
      filter.annotation <- cbind(filter.annotation, "1000gEAS")
    }
    if(row_data[["X1000g2015feb_sas"]] > 0.01){
      filter.annotation <- cbind(filter.annotation, "1000gSAS")
    }
    if(row_data[["X1000g2015feb_afr"]] > 0.01){
      filter.annotation <- cbind(filter.annotation, "1000gAFR")
    }
    if(as.numeric(row_data[["ALT.COUNT"]]) < alt.count & row_data[["cosmic70"]]=="."){
      filter.annotation <- cbind(filter.annotation, paste("ALTCOUNT<", alt.count, sep=""))
    }
    if(row_data[["IS_CANONICAL_TRANSCRIPT"]]=="N"){
      filter.annotation <- cbind(filter.annotation, "NotCanonicalTrans")
    }
    if(is.null(filter.annotation)) {
      # if variant passes all filters, variant is assigned a PASS tag
      filter.annotation <- "PASS"
    }    
    return(paste(filter.annotation, collapse=','))    
  }
  
  # fix 1000g data, where . and numbers are mixed
  variants[variants$X1000g2015feb_eur==".",]$X1000g2015feb_eur <- 0
  variants[variants$X1000g2015feb_amr==".",]$X1000g2015feb_amr <- 0
  variants[variants$X1000g2015feb_eas==".",]$X1000g2015feb_eas <- 0
  variants[variants$X1000g2015feb_sas==".",]$X1000g2015feb_sas <- 0
  variants[variants$X1000g2015feb_afr==".",]$X1000g2015feb_afr <- 0
  
  variants$X1000g2015feb_eur <- as.numeric(variants$X1000g2015feb_eur)
  variants$X1000g2015feb_amr <- as.numeric(variants$X1000g2015feb_amr)
  variants$X1000g2015feb_eas <- as.numeric(variants$X1000g2015feb_eas)
  variants$X1000g2015feb_sas <- as.numeric(variants$X1000g2015feb_sas)
  variants$X1000g2015feb_afr <- as.numeric(variants$X1000g2015feb_afr)
  
  # extract all variants that occur in healthy variants with minimum vaf
  # and filter for these
  healthy.variants <- variants[variants$PLOTNAME %in% healthy.samples & variants$VARIANT_FREQUENCY>=vaf.healthy,]$CANONICAL_TRANSCRIPT_AAChange
  
  #ALT.COUNT should be at larger than 2
  variants$ALT.COUNT <- sapply(variants$AD, function(x) strsplit(as.character(x), ",")[[1]][2])
  
  return(apply(variants, 1, .Annotate))
}
