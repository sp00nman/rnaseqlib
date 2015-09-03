library(ggplot2)
library(reshape)
library(pheatmap)
library(RColorBrewer)
library(scales)


read.clinical.info <- function(clinical_info,
                               header=TRUE){

    patients <- read.table(
        file=clinical_info,
        sep="\t",
        header=header,
        stringsAsFactors=FALSE,
        na.strings=""
    )

    return(patients)
}


read.mutation.info <- function(mutation_info, header=FALSE){

    REF.COUNT <- c()
    ALT.COUNT <- c()

    mutations <- read.table(
        file=mutation_info,
        sep="\t",
        header=header,
        stringsAsFactors=FALSE,
        na.strings=""
    )

    colnames(mutations) <- c("UNIQ_PATIENT_ID",
                             "PATIENT_ID",
                             "ENSGENEID",
                             "AD",
                             "GENESYMBOL")

    ref.alt <- strsplit(mutations$AD, ",")

    for (i in 1:length(ref.alt)){
        REF.COUNT <- cbind(REF.COUNT,
                           ref.alt[[i]][1])
        ALT.COUNT <- cbind(ALT.COUNT,
                           ref.alt[[i]][2])
    }

    mutations$REF.COUNT <- as.numeric(REF.COUNT)
    mutations$ALT.COUNT <- as.numeric(ALT.COUNT)
    mutations$AF <- mutations$ALT.COUNT/(mutations$ALT.COUNT + mutations$REF.COUNT)

    return(mutations)
}


exclude.diagnosis.datasets <- function(dta, exclusion) {

    return(dta[dta$DIAGNOSIS!=exclusion,])
}


create.matrix <- function(dta,column="AF"){

  num_nrow <- length(unique(dta$SAMPLE_ID))
  num_col <- length(unique(na.omit(dta$GENESYMBOL)))

  mut_matrix <- matrix(0, nrow=num_nrow, ncol=num_col)

  rownames(mut_matrix) <- unique(dta$SAMPLE_ID)
  colnames(mut_matrix) <- unique(na.omit(dta$GENESYMBOL))

  for (i in 1:length(dta$SAMPLE_ID)){
    if(is.na(dta$GENESYMBOL[i])){next}
    else {mut_matrix[dta$SAMPLE_ID[i], dta$GENESYMBOL[i]] = dta[[column]][i]}
  }

  return(mut_matrix)
}

count.recurrent.mutations <- function(drivers, mutations, patmut) {
  
  mutation_list <- unique(c(drivers,mutations$GENESYMBOL))
  fill_sum <- vector(length=length(mutation_list))
  
  # TODO: replace the loop with an apply function
  for (i in 1:length(mutation_list) ) {
    gene <- mutation_list[i]
    fill_sum[i] <- sum(patmut[[gene]]!=0, 
                       na.rm=TRUE)
  }
  
  bp.df <- data.frame(mutation_list, 
                      fill_sum)
  colnames(bp.df) <- c("gene", "occurence")
  bp.df$gene <- factor(bp.df$gene, 
                       levels=mutation_list)
  levels.bp <- bp.df[order(bp.df$occurence, 
                           decreasing=FALSE),]$gene
  bp.df$gene <- factor(bp.df$gene, 
                       levels=levels.bp)
  
  # first order by drivers, then all the rest
  driver_occurence <- bp.df[bp.df$gene %in% drivers,]
  mutation_occurence <- bp.df[!bp.df$gene %in% drivers,]
  
  dta_reord <- rbind(driver_occurence[order(driver_occurence$occurence, 
                                            decreasing=TRUE),],
                     mutation_occurence[order(mutation_occurence$occurence, 
                                              decreasing=TRUE),])
  dta_reord$gene <- factor(dta_reord$gene, levels=dta_reord$gene)
  
  return(dta_reord)
}

plot.barplot.occurences <- function(bp.df,revert=TRUE){
  
  if(revert){
    bp.df <- bp.df[dim(bp.df)[1]:1,]
    bp.df$gene <- factor(bp.df$gene, levels=bp.df$gene)
  }
  
  p <- ggplot(data=bp.df, aes(x=gene,y=occurence))
  p = p + geom_bar(colour="white", fill="#4092BD", stat="identity" )
  p = p + scale_x_discrete(expand = c(0, 0)) 
  p = p + scale_y_continuous(expand = c(0, 0))
  p = p + labs(title = "", y = "Number of individuals", x = "")
  p = p + coord_flip()
  p = p + theme_bw()
  p = p + theme(panel.grid.major.x=element_blank(), 
                panel.grid.major.y=element_blank(),
                panel.grid.minor.x=element_blank(), 
                panel.grid.minor.y=element_blank(),
                panel.border = element_blank(), 
                axis.line = element_line(colour="black"),
                axis.text.y =  element_text(face = "italic", 
                                            color = "black"),
                legend.position="none")
  return(p)
}


print_plot <- function(filename, g_obj){
  pdf(file=filename, width=14, height=10, useDingbats=FALSE)
  print(g_obj)
  dev.off()
}

create.patient.mutation.matrix <- function(patmut, bp.df, drivers) {
  
  # order matrix, first order by diagnosis
  # then by number of drivers
  # then rest
  fill_matrix <- matrix(nrow=dim(patmut)[1], ncol=length(bp.df$gene)+1)
  fill_matrix[,1] <- patmut$DIAGNOSIS
  count <- length(bp.df$gene) + 1 
  
  for (i in 2:count) {
    gene <- bp.df$gene[i-1]
    fill_matrix[,i] <- patmut[[as.character(gene)]]
  }
  
  ord.df <- as.data.frame(fill_matrix)
  colname_df <- bp.df$gene
  colnames(ord.df) <- c("DIAGNOSIS", as.character(colname_df))
  
  # Sort by all columns in the data frame, from left to right
  ord_mut <- patmut[do.call(order, c(as.list(ord.df), decreasing=TRUE)), ]
  levels <-ord_mut$SAMPLE_ID
  ord_mut$SAMPLE_ID <- factor(ord_mut$SAMPLE_ID, levels=levels)
  
  #creat dataframe for plotting
  fill_mat <- matrix(ncol=length(bp.df$gene), nrow=dim(ord_mut)[1])
  for (i in 1:length(bp.df$gene) ) {
    gene <- bp.df$gene[i]
    #print(gene)
    fill_mat[,i] <-ord_mut[[as.character(gene)]]
  }
  
  # fill dataframe
  fill_df <- as.data.frame(fill_mat)
  plot_m <- data.frame(ord_mut$SAMPLE_ID, fill_df)
  colnames(plot_m) <- c("SAMPLE_ID", as.character(bp.df$gene))
  plot_m$SAMPLE_ID <- sub("H_", "", plot_m$SAMPLE_ID)
  plot_m$SAMPLE_ID <- factor(plot_m$SAMPLE_ID, levels=plot_m$SAMPLE_ID)
  
  # reshape data
  heat.df<-melt(plot_m)
  colnames(heat.df) <- c("SAMPLE_ID", "gene", "count")
  heat.df$gene <- factor(heat.df$gene, levels=rev(as.character(bp.df$gene)))
  
  return(heat.df)
}

plot.patient.mutation.matrix <- function(heat.df){
  
  mut_mat <- ggplot(heat.df, 
                    aes(x=SAMPLE_ID, 
                        y=gene))
  mut_mat <- mut_mat + geom_tile(aes(fill=count), 
                                 colour="white")
  mut_mat <- mut_mat + scale_fill_gradient2(
    low="#ebebe1ff",
    mid=muted("red"), 
    high=muted("red"), 
    midpoint=1, 
    na.value="#ebebe1ff")
  #mut_mat <- mut_mat + scale_fill_manual(
  # values = c("#F3F2ED","#F0A257"), 
  # na.value="#D1D2D4")
  mut_mat <- mut_mat + theme(
    axis.text.x=element_text(angle = -90, 
                             hjust = 0, 
                             size=7, 
                             color="black"),
    axis.ticks.x=element_blank(), 
    axis.ticks.y=element_blank(),
    axis.text.x=element_blank(), 
    axis.text.y=element_text(face = "italic",color = "black"), 
    axis.title.x=element_text(face="bold"), 
    axis.title.y=element_blank(),
    panel.grid.major.x=element_blank(), 
    panel.grid.major.y=element_blank(),
    panel.grid.minor.x=element_blank(), 
    panel.grid.minor.y=element_blank(), 
    panel.background=element_rect(fill="#ffffff")
  )
  
  return(mut_mat)
  
}

options <- commandArgs(trailingOnly=TRUE)
clinical_info = options[1]
mutation_info = options[2]
output_file = options[3]

#clinical_info = "~/Dropbox/CeMM_PhD_Sync2Dropbox/projects_cemm/prj_rnaseq_varcall_myeloid_panel/R_mutation_patient_plot/patient_matrix.csv"
#mutation_info = "~/Dropbox/CeMM_PhD_Sync2Dropbox/projects_cemm/prj_rnaseq_varcall_myeloid_panel/R_mutation_patient_plot/mut_patient_mutation.table"

# read in data
patients <- read.clinical.info(clinical_info)
mutations <- read.mutation.info(mutation_info)
# merge datasets
merge_dta <- merge(
  patients,
  mutations,
  by.x="SAMPLE_ID",
  by.y="PATIENT_ID",
  all.x=TRUE
)
# exclude healthy [HEALTHY_1] patients for plotting
excluded_patients <- exclude.diagnosis.datasets(
  merge_dta, 
  exclusion=1
)
# create matrix
patients_mut_mat <- create.matrix(excluded_patients)
# convert matrix to dataframe and merge again
dta_patmut <- as.data.frame(patients_mut_mat)
dta_patmut$UNIQ_SAMPLE_ID <- rownames(dta_patmut)
merge_patmut <- merge(
  excluded_patients, 
  dta_patmut, 
  by.x="SAMPLE_ID", 
  by.y="UNIQ_SAMPLE_ID"
)

uniq_merge_patmut <- merge_patmut[!duplicated(merge_patmut[c("SAMPLE_ID")]), ]
# count reccurent mutations
bp.df <- count.recurrent.mutations(
  drivers=c("JAK2_CLIN",
            "CALR_CLIN",
            "MPL_CLIN"),
  mutations,
  patmut=uniq_merge_patmut
)

g_obj <- plot.barplot.occurences(bp.df)
print_plot(output_file, 
           g_obj)

plot_matrix <- create.patient.mutation.matrix(
  patmut = uniq_merge_patmut,
  bp.df = bp.df,
  drivers = drivers
)
j_obj <- plot.patient.mutation.matrix(plot_matrix)
print_plot(output_file, 
           j_obj)




