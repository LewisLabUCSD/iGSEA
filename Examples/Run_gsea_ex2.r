rm(list=ls())
setwd("~/GitHub/Bioinformatics_Tools/")
source('R/Gsea.R')

########################################################################################
##  gsea_ex2
#---------------------------------------------------------------------------------------
#' @description This example is designed to demonstrate how to run gsea analysis on the 
#'              input file (DESeq2 object): '.rda' file 
#' @details All the parameter definitions is in the below section of "Setting Parameter"
#' @return out a list that summarizes the GSEA enrichment analysis
#' 
#' @export exportExcel_fp an excel file of the 'out' can be exported
#' @export exportPDF_fp an figure of the GSEA enrichment analysis
#' --------------------------------------------------------------------------------------------------
#' @author Austin W.T. Chiang
#' @note   v1.4 -- 2017.02.15 Add description of ex1, ex2, and idmapping file
#' @note   v1.3 -- 2017.02.08 Add function of displaying topN pathways and optimze code
#' @note   v1.2 -- 2017.02.06 Add zebrafish2 human id mapping
#' @note   v1.1 -- 2017.01.19 Midified for BEN
#' @note   v1.0 -- 2016.12.15 File created
#' --------------------------------------------------------------------------------------------------
########################################################################################
gsea_ex2 <- function(DEout_obj_fp = NULL, export_path = NULL,
                     IDmappingfile = NULL, convert = FALSE,
                     gene_name = TRUE, GSEA_set = c('c0', 'c2', 'c3'),
                     minGSSize = 30, pAdjustMethod = "none",
                     pvalueCutoff = 0.01, export = "BOTH",
                     ShowTop = TRUE, TopN = 5){
  print(paste0("Here we go .... "))
  ## Input the DESeq result file
  x <- load(DEout_obj_fp)
  myClass <- get(x)
  if(convert == TRUE){
      geneMapTable <- read.csv(IDmappingfile, sep = "\t")
      myClass <- useIDmappingFile(myClass=myClass,
                                  geneMapTable=geneMapTable)
  }
  
  #DEGTable <- subset(myClass, abs(log2FoldChange) >= .58 & (padj < 0.05) )
  DEGTable <- subset(myClass, (padj < 0.05) )
  DEGTable <- DEGTable[with(DEGTable, order(-log2FoldChange)), ]
  if(gene_name==TRUE){
    DEGTable$Gene <- DEGTable$gene_name
  }else{
    DEGTable$Gene <- row.names(DEGTable)
  }
  DEGTable <- DEGTable[,c("Gene","log2FoldChange")]
  
  ## do the GSEA analysis
  exportExcel_fp = paste0(export_path,"gsea_out.xlsx")
  exportPDF_fp = paste0(export_path,"gsea_out.pdf")
  GSEA_res <-  AnnotDEG_byGSEA(DEGTable = DEGTable, GSEA_set = GSEA_set,
                               export = export, ShowTop = ShowTop, TopN = TopN,
                               exportExcel_fp = exportExcel_fp,
                               exportPDF_fp = exportPDF_fp,
                               minGSSize = minGSSize, 
                               pAdjustMethod = pAdjustMethod, 
                               pvalueCutoff = pvalueCutoff)
}

########################################################################################
##  Setting Parameter
########################################################################################
#' @param DEout_obj_fp the RDA file of the DESeq2 output object 
#' @param export_path the path that user want to export the GSEA results
#' @param convert 'TRUE': need to do id converion; 'FALSE': no id conversion needed  
#' @param IDmappingfile the idmapping file show that user's id to human gene name
#' @param gene_name 'TRUE': there is a columan "gene_name" in the input files (i.e. DESeq2 output files)
#' @param GSEA_set = c('c0', 'c2', 'c3')
#'                 co - hallmark gene sets;       c1 - positional gene sets; 
#'                 c2 - curated gene sets;        c3 - motif gene sets;
#'                 c4 - computational gene sets;  c5 - GO gene sets;
#'                 c6 - oncogenic signatures;     c7 - immunologic signatures. 
#'                 (http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C2)
#' @param minGSSize minimum size of genes for a specific GSEA gene set
#' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param pvalueCutoff cutoff of pvalue in the GSEA enrichment analysis 
#' @param export option for exporting GSEA result: "EXCEL", "PDF", or "BOTH"
#' @param ShowTop 'TRUE': the visualization will display only the top and bottom N gene sets; 
#'                'FALSE': all of the significant gsea gene sets will be displayed 
#'                        --> Note that: it is usually results in too many gene sets, which are not easy to see anything    
#' @param TopN this parameter indicate how many gene sets need to be displayed in the plots
########################################################################################

## Input the DESeq result object
DEout_obj_fp <- 'Data/DESeqObj_ex2.rda'
export_path <- "Results/DE_out_test/"
IDmappingfile <- "Data/cho2human_idmapping.txt"
#IDmappingfile <- "Data/mouse2human_idmapping.txt"
#IDmappingfile <- "Data/zebrafish2human_idmapping.txt"
convert = TRUE
gene_name = TRUE 
GSEA_set = c('c0', 'c2', 'c3')
minGSSize = 30
pAdjustMethod = "none"
pvalueCutoff = 0.01
export = "BOTH"
ShowTop = TRUE 
TopN = 5


########################################################################################
##  RUN gsea_ex2
########################################################################################
gsea_ex2(DEout_obj_fp = DEout_obj_fp,    export_path = export_path,     
         IDmappingfile = IDmappingfile,  convert = convert,
         GSEA_set = GSEA_set,            minGSSize = minGSSize,
         pAdjustMethod = pAdjustMethod,  pvalueCutoff = pvalueCutoff,
         export = export,                ShowTop = ShowTop,   
         TopN = TopN)  
  
