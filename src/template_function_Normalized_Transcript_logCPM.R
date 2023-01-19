Normalized_Transcript_logCPM <- function(CCBR_1188_Isoforms_Rawcounts, Normalized_Counts_Isoforms) {
    
    library(tidyverse)
    df <- CCBR_1188_Isoforms_Rawcounts
    df.n <- Normalized_Counts_Isoforms

    genes <- df$GeneName
    names(genes) <- df$transcript_id

    df.n %>% mutate(Genename = genes[Gene]) -> df.n
    df.n %>% select(Genename,c(1:21)) -> df.n
    colnames(df.n)[2] <- "Transcript_ID"
    head(df.n)
    
    return(df.n)
}

print("template_function_Normalized_Transcript_logCPM.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_CCBR_1188_Isoforms_Rawcounts<-readRDS(paste0(rds_output,"/var_CCBR_1188_Isoforms_Rawcounts.rds"))
Input_is_Seurat_count <- 0
for(item in var_CCBR_1188_Isoforms_Rawcounts){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_CCBR_1188_Isoforms_Rawcounts<-as.data.frame(var_CCBR_1188_Isoforms_Rawcounts)}else{var_CCBR_1188_Isoforms_Rawcounts <- var_CCBR_1188_Isoforms_Rawcounts}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Normalized_Counts_Isoforms<-readRDS(paste0(rds_output,"/var_Normalized_Counts_Isoforms.rds"))
Input_is_Seurat_count <- 0
for(item in var_Normalized_Counts_Isoforms){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Normalized_Counts_Isoforms<-as.data.frame(var_Normalized_Counts_Isoforms)}else{var_Normalized_Counts_Isoforms <- var_Normalized_Counts_Isoforms}
invisible(graphics.off())
var_Normalized_Transcript_logCPM<-Normalized_Transcript_logCPM(var_CCBR_1188_Isoforms_Rawcounts,var_Normalized_Counts_Isoforms)
invisible(graphics.off())
saveRDS(var_Normalized_Transcript_logCPM, paste0(rds_output,"/var_Normalized_Transcript_logCPM.rds"))
