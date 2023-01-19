# DEG Analysis [CCBR] (de953b9c-b7f3-49ac-b883-55070d759e5d): v175
DEG_Analysis_Genes <- function(Filtered_Counts_Genes, CCBR_1188_Metadata_Final) {
    
    ## --------- ##
    ## Libraries ##
    ## --------- ##
    
    library(limma)
    library(tidyverse)
    library(edgeR)
    library(stringr)
    library(grid)
    library(gridExtra)
    

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##

    #Basic Parameters:
    counts_matrix <- Filtered_Counts_Genes 
    sample_metadata <- CCBR_1188_Metadata_Final
    gene_names_column="Gene"
    sample_name_column<-"Sample"
    samples_to_include = c("JINB8_CD3_CD28_1","JINB8_CD3_CD28_3","JINB8_CD3_CD28_6","JINB8_CD3_CD28_TSP1_1","JINB8_CD3_CD28_TSP1_3","JINB8_CD3_CD28_TSP1_6","JINB8_TSP1_1","JINB8_TSP1_3","JINB8_TSP1_6","JINB8_UT_1","Jurkat_CD3_CD28_1","Jurkat_CD3_CD28_3","Jurkat_CD3_CD28_6","Jurkat_CD3_CD28_TSP1_1","Jurkat_CD3_CD28_TSP1_3","Jurkat_CD3_CD28_TSP1_6","Jurkat_TSP1_1","Jurkat_TSP1_3","Jurkat_TSP1_6","Jurkat_UT_1")
    contrast_variable_column<-"Cell_Line"
    contrasts<-c("JINB8-JK")
    covariates_columns=c("Cell_Line","Batch")

    #Advanced Parameters:
    input_in_log_counts <- FALSE
    return_mean_and_sd<-FALSE
    return_normalized_counts<-TRUE
    normalization_method<-"quantile"
    
    ##--------------- ##
    ## Error Messages ##
    ## -------------- ##

    if(make.names(colnames(counts_matrix))!=colnames(counts_matrix)){
        print("Error: The following counts matrix column names are not valid:\n")
        print(colnames(counts_matrix)[make.names(colnames(counts_matrix))!=colnames(counts_matrix)])

        print("Likely causes are columns starting with numbers or other special characters eg spaces.")
        stop("Bad column names.")
    }

    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##
    
    samples_to_include <- samples_to_include[samples_to_include != "Gene"]
    samples_to_include <- samples_to_include[samples_to_include != "GeneName"]
    df.m <- counts_matrix[,samples_to_include]
    gene_names <- NULL
    gene_names$GeneID <- counts_matrix[,gene_names_column]
    
    ### This code block does input data validation
    sample_metadata <- sample_metadata[match(colnames(df.m),sample_metadata[,sample_name_column]),]
    sample_metadata <- sample_metadata[rowSums(is.na(sample_metadata)) != ncol(sample_metadata), ]
    df.m <- df.m[,match(sample_metadata[,sample_name_column],colnames(df.m))]
    
    #Create DGEList object from counts
    if(input_in_log_counts == TRUE){
        x <- DGEList(counts=2^df.m, genes=gene_names)
    } else {
        x <- DGEList(counts=df.m, genes=gene_names) 
    }
    
    #Put covariates in order 
    covariates_columns=covariates_columns[order(covariates_columns!=contrast_variable_column)]

    for(ocv in covariates_columns){
        sample_metadata[,ocv]=gsub(" ","_",sample_metadata[,ocv])
    }

    contrasts=gsub(" ","_",contrasts)

    dm.formula <- as.formula(paste("~0 +", paste(covariates_columns, sep="+", collapse="+")))
    design=model.matrix(dm.formula, sample_metadata)
    colnames(design) <- str_replace_all(colnames(design), contrast_variable_column, "")
    
    if (normalization_method %in% c("TMM","TMMwzp","RLE","upperquartile")){
        x <- calcNormFactors(x, method = normalization_method) 
        rownames(x) <- x$genes$GeneID
        v <- voom(x,design=design,normalize="none")
    } else {
        v <- voom(x,design=design,normalize=normalization_method,plot=TRUE,save.plot = TRUE)
    }
    
    rownames(v$E) <- v$genes$GeneID
    as.data.frame(v$E) %>% rownames_to_column("Gene") -> df.voom
    fit <- lmFit(v, design)

    #Print Mean-variance Plot
    sx <- v$voom.xy$x
    sy <- v$voom.xy$y
    xyplot <- as.data.frame(cbind(sx,sy))
    voomline <- as.data.frame(cbind(x=v$voom.line$x,y=v$voom.line$y))
    
    g <- ggplot() +
        geom_point(data=xyplot, aes(x=sx,y=sy),size=1) +
        theme_bw() +
        geom_smooth(data=voomline, aes(x=x,y=y),color = "red") +
        ggtitle("voom: Mean-variance trend") +
        xlab(v$voom.xy$xlab) + ylab(v$voom.xy$ylab) + 
        theme(axis.title=element_text(size=12),
        plot.title = element_text(size = 14, face = "bold",hjust = 0.5))

    print(g)

    #Run Contrasts
    cm <- makeContrasts(contrasts = contrasts, levels=design)
    
    fit2 <- contrasts.fit(fit, cm)
    fit2 <- eBayes(fit2)
    logFC = fit2$coefficients
    colnames(logFC)=paste(colnames(logFC),"logFC",sep="_")
    tstat = fit2$t
    colnames(tstat)=paste(colnames(tstat),"tstat",sep="_")
    FC = 2^fit2$coefficients
    FC = ifelse(FC<1,-1/FC,FC)
    colnames(FC)=paste(colnames(FC),"FC",sep="_")
    pvalall=fit2$p.value
    colnames(pvalall)=paste(colnames(pvalall),"pval",sep="_")
    pvaladjall=apply(pvalall,2,function(x) p.adjust(x,"BH"))
    colnames(pvaladjall)=paste(colnames(fit2$coefficients),"adjpval",sep="_")

    
    if(return_mean_and_sd == TRUE){
        tve <- t(v$E)        
        mean.df <- as.data.frame(tve) %>% rownames_to_column("Sample") %>% mutate(group=sample_metadata[sample_metadata[,sample_name_column]==Sample,contrast_variable_column]) %>% group_by(group) %>% summarise_all(funs(mean)) %>% as.data.frame()
        mean.df[,-c(1,2)] %>% as.matrix() %>% t() -> mean
        colnames(mean) <- mean.df[,1]
        colnames(mean)=paste(colnames(mean),"mean", sep="_")
        colnames(mean) = gsub("\\.", "_", colnames(mean))
        
        sd.df <- as.data.frame(tve) %>% rownames_to_column("Sample") %>% mutate(group=sample_metadata[sample_metadata[,sample_name_column]==Sample,contrast_variable_column]) %>% group_by(group) %>% summarise_all(funs(sd)) %>% as.data.frame()
        sd.df[,-c(1,2)] %>% as.matrix() %>% t() -> sd
        colnames(sd) <- sd.df[,1]
        colnames(sd)=paste(colnames(sd), "sd",sep="_")
        colnames(sd) = gsub("\\.", "_", colnames(sd))
    finalres=as.data.frame(cbind(mean, sd,  FC, logFC, tstat, pvalall, pvaladjall)) 
    } else {
        finalres=as.data.frame(cbind(FC, logFC, tstat, pvalall, pvaladjall))
    }

    if(return_normalized_counts == TRUE){
        finalres = as.data.frame(cbind(finalres, v$E))
    }

    finalres %>% rownames_to_column("Gene") -> finalres
    print(paste0("Total number of genes included: ", nrow(finalres)))
    
    getgenelists <- function(FClimit,pvallimit,pval){
        upreggenes <- list()
        downreggenes <- list()
        for(i in 1:length(contrasts)){
            if(pval == "pval"){
            finalres %>% dplyr::filter(.data[[colnames(FC)[i]]] > FClimit & .data[[colnames(pvalall)[i]]] < pvallimit) %>% pull(Gene) %>% length() -> upreggenes[[i]] 
            finalres %>% dplyr::filter(.data[[colnames(FC)[i]]] < -FClimit & .data[[colnames(pvalall)[i]]] < pvallimit) %>% pull(Gene) %>% length() -> downreggenes[[i]]        
        } else {
            finalres %>% dplyr::filter(.data[[colnames(FC)[i]]] > FClimit & .data[[colnames(pvaladjall)[i]]] < pvallimit) %>% pull(Gene) %>% length() -> upreggenes[[i]] 
            finalres %>% dplyr::filter(.data[[colnames(FC)[i]]] < -FClimit & .data[[colnames(pvaladjall)[i]]] < pvallimit) %>% pull(Gene) %>% length() -> downreggenes[[i]] 
        }
        }
        names(upreggenes) <- contrasts
        names(downreggenes) <- contrasts
        allreggenes <- rbind(unlist(upreggenes),unlist(downreggenes))
        rownames(allreggenes) <- c(paste0("upreg>",FClimit, ", ",pval,"<",pvallimit),paste0("downreg<-",FClimit, ", ",pval,"<",pvallimit))
        return(allreggenes)
    }

    FCpval1 <- getgenelists(FClimit = 1.2, pvallimit = 0.05,"pval")
    FCpval2 <- getgenelists(FClimit = 1.2, pvallimit = 0.01,"pval")
    FCadjpval1 <- getgenelists(FClimit = 1.2, pvallimit = 0.05,"adjpval")
    FCadjpval2 <- getgenelists(FClimit = 1.2, pvallimit = 0.01,"adjpval")

    pvaltab <- rbind(FCpval1,FCpval2,FCadjpval1,FCadjpval2)
    table <- tableGrob(pvaltab)
    grid.newpage()
    grid.draw(table)

    grid.arrange(g,table, ncol = 1)

    call_me_alias<-colnames(finalres)
    colnames(finalres)<-gsub("\\(|\\)","", call_me_alias)
    df.final<-(finalres)
   
    return(df.final) 
}

print("template_function_DEG_Analysis_Genes.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Filtered_Counts_Genes<-readRDS(paste0(rds_output,"/var_Filtered_Counts_Genes.rds"))
Input_is_Seurat_count <- 0
for(item in var_Filtered_Counts_Genes){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Filtered_Counts_Genes<-as.data.frame(var_Filtered_Counts_Genes)}else{var_Filtered_Counts_Genes <- var_Filtered_Counts_Genes}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_CCBR_1188_Metadata_Final<-readRDS(paste0(rds_output,"/var_CCBR_1188_Metadata_Final.rds"))
Input_is_Seurat_count <- 0
for(item in var_CCBR_1188_Metadata_Final){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_CCBR_1188_Metadata_Final<-as.data.frame(var_CCBR_1188_Metadata_Final)}else{var_CCBR_1188_Metadata_Final <- var_CCBR_1188_Metadata_Final}
invisible(graphics.off())
var_DEG_Analysis_Genes<-DEG_Analysis_Genes(var_Filtered_Counts_Genes,var_CCBR_1188_Metadata_Final)
invisible(graphics.off())
saveRDS(var_DEG_Analysis_Genes, paste0(rds_output,"/var_DEG_Analysis_Genes.rds"))
