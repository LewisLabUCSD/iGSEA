library(ggplot2)
library(GSEABase)
library(stringr)
library(clusterProfiler)
#library(WriteXLS)
library(gdata)
library(data.table)
library(DOSE)
library(openxlsx)
library(biomaRt)

#' Annotate_GSEA:
#'    This function is used to annotate the input gene list with GSEA enrichment analysis.
#' --------------------------------------------------------------------------------------------------
#' @param DEGTable a table of DEGs with Human gene names 
#' @param pvalueCutoff cutoff of pvalue in the GSEA enrichment analysis 
#' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param minGSSize minimum size of genes for a specific GSEA gene set
#' @param ShowTop whether only show the top&bottom pathways
#' @param TopN how many top&bottom 'N' pathways will be display
#' @param qvalueCutoff minimum size of genes for GSEA analysis
#' @param universe universe background genes
#' @param export option for exporting GSEA result: "EXCEL", "PDF", or "BOTH"
#' @param GSEA_set a list of interested GSEA sets:
#'                 co - hallmark gene sets;       c1 - positional gene sets; 
#'                 c2 - curated gene sets;        c3 - motif gene sets;
#'                 c4 - computational gene sets;  c5 - GO gene sets;
#'                 c6 - oncogenic signatures;     c7 - immunologic signatures. 
#'                 (http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C2)
#'
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
#' @note   v1.0 -- 2016.12.15
#' --------------------------------------------------------------------------------------------------
#'
AnnotDEG_byGSEA <- function(DEGTable = NULL, 
                            minGSSize = 30, ShowTop = TRUE, TopN = 10,
                            pAdjustMethod = "none", pvalueCutoff = 0.01,
                            export = "BOTH", exportExcel_fp = "GSEA_out.xlsx", exportPDF_fp = "GSEA_out.pdf",
                            GSEA_set = NULL){
  geneList <- DEGTable[,"log2FoldChange"]
  names(geneList) <- toupper(DEGTable[,"Gene"])

  geneList <- na.omit(geneList)
  geneList1 <- sort(geneList,decreasing = TRUE)  
  
  ## SKIP the analysis if the gene list is too short (less than 10 genes)
  out <- list()
  if (length(geneList1) <= 30)
    return(out)
  
  ## do the GSEA enrichment analysis
  egmt <- list()
  for (c in GSEA_set) {
    gmtfile_name <-
      paste0("msigdb_v5.1_GMTs/", c, ".all.v5.1.symbols.gmt")
    gmt <- read.gmt(gmtfile_name)
    egmt[[c]]     <- tryCatch(  GSEA(geneList=geneList1, 
                                     TERM2GENE = gmt,
                                     pvalueCutoff = pvalueCutoff,
                                     pAdjustMethod = pAdjustMethod, 
                                     minGSSize = minGSSize ), error=function(e) NULL)
  }
  
  ## summarize the GSEA enrichment analysis for output
  for (c in names(egmt)) {
    tmp = NULL
    if (!is.na(egmt[[c]])) {
      tmp <- summary(egmt[[c]])
      if((dim(tmp)[1] > 0) && (tmp$setSize >= minGSSize)){
        out[[c]] <- tmp
      }
    }
  }
  
  if(length(out) == 0)
    return(out)
    
  ## summarize the GSEA enrichment analysis for visualize
  count = 1
  for (c in names(out)){
    tmp <- subset(out[[c]], select = c(NES))
    switch(c,
           c0 = {sigName = "Hallmark"},
           c1 = {sigName = "Positional"},
           c2 = {sigName = "Pathway"},
           c3 = {sigName = "TFBS"},
           c4 = {sigName = "Comp. gene set"},
           c5 = {sigName = "Gene Ontology"},
           c6 = {sigName = "Oncogenic Signatures"},
           c7 = {sigName = "Immunological Signatures"})
    tmp$Signature <- sigName
    if(count == 1){
      y <- tmp
      count = 2
    } else {
      y <- rbind(y, tmp)
    }
  }
  colnames(y) <- c("DataSet","Signature")
  wth_factor = dim(y)[1]/300
  
  ## export the GSEA enrichment analysis results
  switch(export,
         EXCEL = ExportExcel(out, exportExcel_fp),
         PDF = ExportPDF(y, exportExcel_fp,wth_factor, ShowTop, TopN),
         BOTH = {ExportExcel(out, exportExcel_fp)
                 ExportPDF(y, exportPDF_fp,wth_factor, ShowTop, TopN)}  )

  return(out)
}

#' ExportExcel
#'
#' @param out 
#' @param exportExcel 
#'
ExportExcel <- function(out, exportExcel){
  ## Write the Excel file
  options(warn=-1)
  if(length(out) > 0){
    #WriteXLS(out, exportExcel, AdjWidth = FALSE,BoldHeaderRow = TRUE, row.names = F)
    write.xlsx(out, exportExcel, AdjWidth = FALSE,BoldHeaderRow = TRUE, row.names = F)
  }
  options(warn=0)
}


#' read.gmt : This function loads the GSEA gene sets for analyzing
#'
#' @param gmtfile the data sets downloaded from GSEA database 
#'                (http://software.broadinstitute.org/gsea/index.jsp)
#'
#' @return ont2gene Term_vs_Gene 1_to_1 mapping
#'
read.gmt <- function(gmtfile) {
  gmt <- getGmt(con=gmtfile)
  ont2gene <- geneIds(gmt) %>% stack
  ont2gene <- ont2gene[, c("ind", "values")]
  colnames(ont2gene) <- c("ont", "gene")
  
  return(ont2gene)
}


#' get_GSEAPlot
#'
#' @param y 
#' @param fp 
#' @param factor 
#' @param title 
#'
#' @return
#' @export
#'
#' @examples
ExportPDF <- function(y,fp,wth_factor,ShowTop, TopN){
  y.dt <- as.data.table(y)
  y.dt <- y.dt[,names:= rownames(y)]
  melt.1 <- melt(y.dt, id.vars = 'names', measure.vars = c('DataSet'))
  melt.1 <- melt.1[(names != ""),] 
  melt.2 <- melt(y.dt, id.vars = 'names', measure.vars = c('Signature'))
  melt.2 <- melt.2[(names != ""),] 
  y.m  <- merge(melt.1, melt.2, by = 'names')
  colnames(y.m) <- c("names","DataSet","NES","variable.y","Signature")
  y.m$NES <- as.numeric(y.m$NES)
  y.m$Signature <- factor(y.m$Signature,
                          levels = c("Hallmark",
                                     "Positional gene sets",
                                     "Pathway",
                                     "TFBS",
                                     "Comp. gene set",
                                     "Gene Ontology",
                                     "Oncogenic Signatures",
                                     "Immunological Signatures"))
  y.m$names = str_wrap( gsub('_',' ',y.m$names),width=50 )
  y.m$names<-reorder(y.m$names, y.m$NES)

  ## Only plot top/bottom N pathways
  if((ShowTop == TRUE) && (TopN > 0)){
    y.f <- NULL
    for(n in unique(y.m$Signature)){
      b <-subset(y.m, Signature == n  )
      if(dim(b)[1] >= 2*TopN){
        h <-head(b[order(NES),],TopN)
        t <-tail(b[order(NES),],TopN)
        y.f <- rbind(y.f,h,t)
      } else {
        y.f <- rbind(y.f,b)
      }
    }
    y.m <- y.f
  }
  
  #pdf(file=fp, width=15, height=15*wth_factor,onefile = TRUE)
  pdf(file=fp, width=15, height=15,onefile = TRUE)
  p <- ggplot(y.m, aes(x=names, y=NES)) +
       geom_bar(stat = 'identity', aes(fill = NES))  +
       scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limits=c(-3.5,3.5)) +
       labs(list(x = "Molecular Signature", y = "Normalized Enrichment Score (NES)")) +
       theme(text = element_text(size=16)) +
       theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
       ylim(-3.5, 3.5) +
       coord_flip() +
       facet_grid(Signature~., scales = "free", space = "free") 
  
  print(p)
  dev.off()

  return(p)
}

#' useBioMarT
#'
#' @param myClass 
#' @param dataset 
#' @param Original_attribute original gene/transcript ID
#' @param Convert_attribute the ids that user want to convert (e.g., gene name)
#'
#' @return myClass with the converted ids attached

useBioMarT <- function(myClass = NULL, dataset = NULL, 
                       Original_attribute = NULL,
                       Convert_attribute = NULL){
  # Get an archived version of ensembl
  ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset = dataset, host="www.ensembl.org")
  
  # example list of ENSG ids
  ensemblID = rownames(myClass)
  
  # use this list to get corresponding Gene Symbols
  attributes = listAttributes(ensembl)
  results = getBM(attributes = c(Convert_attribute,Original_attribute), 
                  filters = Original_attribute, values = ensemblID, mart = ensembl)
  myClass$gene_name <- results$mgi_symbol[match(rownames(myClass),results$ensembl_gene_id)]
  return(myClass)
}

#' useBioMarT
#'
#' @param myClass 
#' @param geneMapTable the table for id conversion 
#'
#' @return myClass with the converted ids attached

useIDmappingFile <- function(myClass = NULL, geneMapTable = NULL){
  geneList   <- row.names(myClass)
  gene_name <- geneMapTable$Human_genename[match(geneList, geneMapTable$ID)]
  myClass$gene_name <- toupper(gene_name)
  myClass <- na.omit(myClass)
  return(myClass)
}

