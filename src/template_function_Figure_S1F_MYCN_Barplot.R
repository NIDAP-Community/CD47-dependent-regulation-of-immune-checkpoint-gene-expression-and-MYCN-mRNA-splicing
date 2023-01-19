Figure_S1F_MYCN_Barplot <- function(Normalized_Transcript_logCPM) {

 library(reshape2)
 library()
 df <- Normalized_Transcript_logCPM
 gene <- "MYCN"
 df %>% filter(Genename == gene) -> df.select

 
 df.melt <- melt(df.select)

    df.melt$variable <- factor(df.melt$variable,levels = c("Jurkat_UT_1","Jurkat_CD3_CD28_1","Jurkat_CD3_CD28_3","Jurkat_CD3_CD28_6","Jurkat_CD3_CD28_TSP1_1","Jurkat_CD3_CD28_TSP1_3","Jurkat_CD3_CD28_TSP1_6","Jurkat_TSP1_1","Jurkat_TSP1_3","Jurkat_TSP1_6","JINB8_UT_1","JINB8_CD3_CD28_1","JINB8_CD3_CD28_3","JINB8_CD3_CD28_6","JINB8_CD3_CD28_TSP1_1","JINB8_CD3_CD28_TSP1_3","JINB8_CD3_CD28_TSP1_6","JINB8_TSP1_1","JINB8_TSP1_3","JINB8_TSP1_6"))          
    df.melt$value <- 2^df.melt$value

 df.melt %>% ggplot(aes(x=variable,y=value,fill=Transcript_ID)) +
            geom_bar(stat="identity", color="black", position=position_dodge()) +
            ylab("Normalized Counts per Million") +
            ggtitle(gene) +
            theme_classic() +
            theme(axis.title.y = element_text(size=22,margin = margin(r = 70)), axis.title.x = element_text(size=22,margin = margin(t = 30)),
                    axis.text.y = element_text(size = 15, colour="black"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour="black", size = 15), 
                    legend.key.size = unit(1, 'cm'), legend.title = element_text(size=20,margin = margin(l = 40)), legend.text=element_text(colour="black", size=15), 
                    legend.position = c(0,1), legend.justification = c(0, 1),legend.box.margin = margin(5, l = 5, unit = "mm"), plot.title = element_text(colour="black", size = 30, hjust = 0.5)) -> gplot

library(png)
 imageWidth = 3000
 imageHeight = 3000
 dpi = 300
 
png(
   filename="MYCN.png",
   width=imageWidth,
   height=imageHeight,
   units="px",
   pointsize=4,
   bg="white",
   res=dpi,
   type="cairo")
 
print(gplot)

dev.off()

# auto removed: return(NULL)

}

print("template_function_Figure_S1F_MYCN_Barplot.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_Normalized_Transcript_logCPM<-readRDS(paste0(rds_output,"/var_Normalized_Transcript_logCPM.rds"))
Input_is_Seurat_count <- 0
for(item in var_Normalized_Transcript_logCPM){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_Normalized_Transcript_logCPM<-as.data.frame(var_Normalized_Transcript_logCPM)}else{var_Normalized_Transcript_logCPM <- var_Normalized_Transcript_logCPM}
invisible(graphics.off())
var_Figure_S1F_MYCN_Barplot<-Figure_S1F_MYCN_Barplot(var_Normalized_Transcript_logCPM)
invisible(graphics.off())
saveRDS(var_Figure_S1F_MYCN_Barplot, paste0(rds_output,"/var_Figure_S1F_MYCN_Barplot.rds"))
