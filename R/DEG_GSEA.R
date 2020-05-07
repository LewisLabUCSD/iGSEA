##DEG analysis
source('DEG_include.R')
library(data.table)
library(xlsx)
library(foreach)
library(doParallel)

## Innate Immune List
library(readxl)
InnateImmuneList <- sapply(read_excel("~/GoogleDrive/FDA_VirusModel/InnateImmuneList.xlsx"), toupper)

## for FC by a different model

pvadj_cutoff <- .05
# pvadj_cutoff <- 1
# lfc_cutoff <- 1.5
lfc_cutoff <- .58
# lfc_cutoff <- 0

id_conversion <-
  fread(
    'cho_human_mouse_id_mapping.txt',
    sep = '\t',
    header = T,
    na.strings = '-'
  )
setkey(id_conversion, cho_id)

library(grid)
library(fgsea)

# convert_to = 'mouse_id'
convert_to = 'human_id'
fgseaRes_comparisons_runs = list()

# cl <- makeCluster(6)
# registerDoParallel(cl)
# foreach(run =1:12, .packages = c('data.table', 'fgsea', 'gridExtra', 'ggplot2', 'grid')) %dopar% {

load("FDA_virus_study/FDA_DE.rda") #DiffExp

##DEG exploration
vsv <- merge(as.data.table(t(DiffExp$VSV$FoldChange), keep.rownames = T), as.data.table(t(DiffExp$VSV$Sigp), keep.rownames = T), by = 'rn', suffixes = c('.lfc', '.p'))
vsv[, mouse_id:=as.integer(rn)]
vsv[, rn:=NULL]
write.csv(vsv, file = 'vsv.csv', row.names = F)
vsv.converted <- merge(vsv, id_conversion[!duplicated(id_conversion[, c('human_id', 'mouse_id'), with = F])], by = 'mouse_id')
write.csv(vsv.converted, file = 'vsv_converted.csv', row.names = F)

for (virus in names(DiffExp)) {
  res <-
    as.data.table((t(DiffExp[[virus]]$FoldChange)), keep.rownames = T)
  res.Sigp <-
    as.data.table((t(DiffExp[[virus]]$Sigp)), keep.rownames = T)
  models <- names(res)[-1]
  names(res)[1] <- 'mouse_id'
  names(res.Sigp)[1] <- 'mouse_id'
  res[, mouse_id := as.integer(mouse_id)]
  res.Sigp[, mouse_id := as.integer(mouse_id)]
  
  res <-
    merge(res, id_conversion[!duplicated(id_conversion[, c('human_id', 'mouse_id'), with = F])], by = 'mouse_id', all.x = T)[!is.na(human_symbol)][!duplicated(human_symbol)] ## id_conversion
  
  res.Sigp <-
    merge(res.Sigp, id_conversion[!duplicated(id_conversion[, c('human_id', 'mouse_id'), with = F])], by = 'mouse_id', all.x = T)[!is.na(human_symbol)][!duplicated(human_symbol)] ## id_conversion
  merge(res, res.Sigp, by = 'human_symbol', suffixes = c('.lfc', '.p'))
  fgseaRes_comparisons = list()
  
  for (comparison in models) {
    DEGs <- res[!is.na(res[[comparison]])&res.Sigp[[comparison]]<pvadj_cutoff&abs(res[[comparison]])>lfc_cutoff]$human_symbol
    comp_df <- res[!is.na(res[[comparison]])]
    sorted_df <- comp_df[order(comp_df[[comparison]])]
    sorted_expr <- setNames(sorted_df[[comparison]],
                            as.character(sorted_df$human_symbol))
    
    fgseaRes <- list()
    for (i in c(0, 2:6)) {
      ## looping through all gmt files    
      gmtfile <- paste0("msigdb_v5.2_GMTs/c", i, ".all.v5.2.symbols.gmt.txt")
      pathways <- gmtPathways(gmt.file = gmtfile)
      fgseaRes_ <- customfGSEA(
        ## customized for more detailed output summary
        pathways = pathways,
        stats = sorted_expr,
        minSize = 2,
        maxSize = 500,
        nperm = 50000,
        nproc = 4,
        DEGs = DEGs,
        Immu = InnateImmuneList
      )
      # fgseaRes_ <- fgsea(pathways = pathways,
      #                                 stats = sorted_expr,
      #                                 minSize=2, ## restore to 1 after debugging finishes
      #                                 maxSize=500,
      #                                 nperm=50000)
      if (nrow(fgseaRes_) == 0)
        next
      
      topPathwaysUp <-
        fgseaRes_[ES > 0,][head(order(pval), n = 10), pathway]
      topPathwaysDown <-
        fgseaRes_[ES < 0,][head(order(pval), n = 10), pathway]
      topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
      # browser()
      gsea_table <- plotGseaTable2(
        pathways = pathways[topPathways],
        stats = sorted_expr,
        fgseaRes = fgseaRes_,
        title = paste(list(
          paste0(virus),
          comparison,
          gmtfile
        ), collapse = '\n'),
        gseaParam = 0.5,
        colwidths = c(6.4, 
                      8, 0.9, 1.2, 1.2,6, 6)
      )
      fn_gsea <-
        paste0(virus,'/', comparison, paste0('_gmt', i))
      # ggsave(paste0('output/gsea_tables_new/', fn_gsea, '.pdf'), gsea_table )
      ggsave(paste0(
        'output/gsea_tables',
        ifelse(pvadj_cutoff >= 1, '/', '_padjcutoff/'),
        fn_gsea,
        '.pdf'
      ),
      gsea_table, scale =2, width = 12, height = 8)
      # ggsave(paste0('output/gsea_tables_padjcutoff/', fn_gsea, '.pdf'), gsea_table )
      fgseaRes_[, totalHits2 := lapply(totalHits, paste0, collapse = ' '), by = pathway]
      fgseaRes_[, leadingEdge2 := lapply(leadingEdge, paste0, collapse = ' '), by = pathway]
      fgseaRes_[, DEGHits2 := lapply(DEGHits, paste0, collapse = ' '), by = pathway]
      fgseaRes_[, ImmuHits2 := lapply(ImmuHits, paste0, collapse = ' '), by = pathway]
      fgseaRes_[, totalHits := NULL]
      fgseaRes_[, leadingEdge := NULL]
      fgseaRes_[, DEGHits := NULL]
      fgseaRes_[, ImmuHits := NULL]
      
      write.csv(
        as.data.table(fgseaRes_),
        file = paste0(
          'output/gsea_tables',
          ifelse(pvadj_cutoff >= 1, '/', '_padjcutoff/'),
          fn_gsea,
          '.csv'
        )
      )
      fgseaRes[[paste0('c', i)]] <-
        list(
          topPathwaysUp = topPathwaysUp,
          topPathwaysDown = topPathwaysDown,
          plot_gsea_table = gsea_table,
          gsea_result = fgseaRes_
        )
    }
    fgseaRes_comparisons[[comparison]] = fgseaRes
  }
  fgseaRes_comparisons_runs[[virus]] = fgseaRes_comparisons
}
save(fgseaRes_comparisons_runs, file = 'fgseaRes_comparisons_runs_1202.Rda')
