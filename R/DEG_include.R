library(ggplot2)
library(fastmatch)
library(gridExtra)
library(BiocParallel)


get_converted_id <- function(cho_id_,
                             lfc_,
                             padj_,
                             id_conversion,
                             convert_to = 'human_id',
                             padj_cutoff = 0.05) {
  cho_id_ <- as.integer(cho_id_)
  converted_ids <-
    unique(na.omit(id_conversion[.(cho_id_)][[convert_to]]))
  if (is.na(padj_)) {
    print ('null padj')
    return(NULL)
  }
  else if ((length(converted_ids) == 0) | (padj_ > padj_cutoff)) {
    return(NULL)
  }
  # print(convert_to)
  return (cbind(
    'cho_id' = rep(cho_id_, length(converted_ids)),
    'log2FoldChange' = rep(lfc_, length(converted_ids)),
    convert_to = converted_ids
  ))
}


customfGSEA <-
  function (pathways,
            stats,
            nperm,
            DEGs = NULL,
            Immu = NULL,
            minSize = 1,
            maxSize = Inf,
            nproc = 0,
            gseaParam = 1,
            BPPARAM = bpparam())
  {
    if (nproc != 0) {
      BPPARAM <- MulticoreParam(workers = nproc)
    }
    minSize <- max(minSize, 1)
    stats <- sort(stats, decreasing = TRUE)
    stats <- abs(stats) ^ gseaParam
    pathwaysFiltered <- lapply(pathways, function(p) {
      as.vector(na.omit(fmatch(p, names(stats))))
    })
    pathwaysSizes <- sapply(pathwaysFiltered, length)
    toKeep <-which(minSize <= pathwaysSizes & pathwaysSizes <= maxSize & pathwaysSizes< length(stats))
    m <- length(toKeep)
    if (m == 0) {
      return(
        data.table(
          pathway = character(),
          pval = numeric(),
          padj = numeric(),
          ES = numeric(),
          NES = numeric(),
          nMoreExtreme = numeric(),
          size = integer(),
          leadingEdge = list()
        )
      )
    }
    pathwaysFiltered <- pathwaysFiltered[toKeep]
    pathwaysSizes <- pathwaysSizes[toKeep]
    K <- max(pathwaysSizes)
    npermActual <- nperm
    gseaStatRes <-
      do.call(
        rbind,
        lapply(
          pathwaysFiltered,
          calcGseaStat,
          stats = stats,
          returnLeadingEdge = TRUE
        )
      )
    leadingEdges <- mapply("[", list(names(stats)), gseaStatRes[,
                                                                "leadingEdge"], SIMPLIFY = FALSE)
    totalHits <- mapply("[", list(names(stats)), pathwaysFiltered, SIMPLIFY = FALSE)
    DEGHits <- sapply(totalHits,function(x, table) x[match(x, table, nomatch = 0) > 0], DEGs) # the %in% function
    ImmuHits <- sapply(totalHits,function(x, table) x[match(x, table, nomatch = 0) > 0], Immu) # the %in% function
    
    pathwayScores <- unlist(gseaStatRes[, "res"])
    granularity <- 1000
    permPerProc <- rep(granularity, floor(npermActual / granularity))
    if (npermActual - sum(permPerProc) > 0) {
      permPerProc <- c(permPerProc, npermActual - sum(permPerProc))
    }
    universe <- seq_along(stats)
    counts <- bplapply(permPerProc, function(nperm1) {
      leEs <- rep(0, m)
      geEs <- rep(0, m)
      leZero <- rep(0, m)
      geZero <- rep(0, m)
      leZeroSum <- rep(0, m)
      geZeroSum <- rep(0, m)
      for (i in seq_len(nperm1)) {
        randSample <- sample.int(length(universe), K)
        if (m == 1) {
          randEsP <- calcGseaStat(
            stats = stats,
            selectedStats = randSample,
            gseaParam = 1
          )
        }
        else {
          randEs <- fgsea::calcGseaStatCumulative(
            stats = stats,
            selectedStats = randSample,
            gseaParam = 1
          )
          randEsP <- randEs[pathwaysSizes]
        }
        leEs <- leEs + (randEsP <= pathwayScores)
        geEs <- geEs + (randEsP >= pathwayScores)
        leZero <- leZero + (randEsP <= 0)
        geZero <- geZero + (randEsP >= 0)
        leZeroSum <- leZeroSum + pmin(randEsP, 0)
        geZeroSum <- geZeroSum + pmax(randEsP, 0)
      }
      data.table(
        pathway = seq_len(m),
        leEs = leEs,
        geEs = geEs,
        leZero = leZero,
        geZero = geZero,
        leZeroSum = leZeroSum,
        geZeroSum = geZeroSum
      )
    }, BPPARAM = BPPARAM)
    counts <- rbindlist(counts)
    leEs = leZero = geEs = geZero = leZeroSum = geZeroSum = NULL
    pathway = padj = pval = ES = NES = geZeroMean = leZeroMean = NULL
    nMoreExtreme = nGeEs = nLeEs = size = NULL
    leadingEdge = NULL
    . = "damn notes"
    pvals <-
      counts[, list(
        pval = min((1 + sum(leEs)) / (1 + sum(leZero)),
                   (1 + sum(geEs)) / (1 + sum(geZero))),
        leZeroMean = sum(leZeroSum) / sum(leZero),
        geZeroMean = sum(geZeroSum) / sum(geZero),
        nLeEs = sum(leEs),
        nGeEs = sum(geEs)
      ), by = .(pathway)]
    pvals[, `:=`(padj, p.adjust(pval, method = "BH"))]
    pvals[, `:=`(ES, pathwayScores[pathway])]
    pvals[, `:=`(NES, ES / ifelse(ES > 0, geZeroMean, abs(leZeroMean)))]
    pvals[, `:=`(leZeroMean, NULL)]
    pvals[, `:=`(geZeroMean, NULL)]
    pvals[, `:=`(nMoreExtreme, ifelse(ES > 0, nGeEs, nLeEs))]
    pvals[, `:=`(nLeEs, NULL)]
    pvals[, `:=`(nGeEs, NULL)]
    pvals[, `:=`(size, pathwaysSizes[pathway])]
    pvals[, `:=`(pathway, names(pathwaysFiltered)[pathway])]
    pvals[, `:=`(leadingEdge, list(leadingEdges))]
    pvals[, `:=`(totalHits, list(totalHits))]
    pvals[, `:=`(DEGHits, list(DEGHits))]
    pvals[, `:=`(ImmuHits, list(ImmuHits))]
    
    pvals <- pvals[]
    pvals
  }




plotGseaTable2 <- function (pathways, stats, fgseaRes, title, gseaParam = 1, colwidths = c(1.4, 
                                                                                           3.2, 0.9, 1.2, 1.2, 5, 5)) 
{
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathways.DEG <- setNames(lapply(fgseaRes[match(names(pathways), pathway), DEGHits], function(p) {
    unname(as.vector(na.omit(match(p, names(statsAdj)))))
  }), names(pathways))
  pathways.Immu <- setNames(lapply(fgseaRes[match(names(pathways), pathway), ImmuHits], function(p) {
    unname(as.vector(na.omit(match(p, names(statsAdj)))))
  }), names(pathways))
  pathways <- lapply(pathways, function(p) {
    unname(as.vector(na.omit(match(p, names(statsAdj)))))
  })
  ps <- lapply(names(pathways), function(pn) {
    p <- pathways[[pn]]
    p.DEG <- pathways.DEG[[pn]]
    p.Immu <- pathways.Immu[[pn]]
    annotation <- fgseaRes[match(pn, fgseaRes$pathway),]
    y.start <- ifelse(p%in%p.DEG|p%in%p.Immu, 0, 0)
    # y.start <- ifelse(p%in%p.DEG|p%in%p.Immu, -statsAdj[p], 0)
    colors = ifelse(p%in%p.DEG&p%in%p.Immu, 'green', ## both: green
                    ifelse(p%in%p.DEG, 'blue', ## deg: blue
                        ifelse(p%in%p.Immu, 'red', '#6c6c6d'))) ## immu: red
    list(textGrob(pn, just = "right", x = unit(0.95, "npc")), 
         ggplot() + geom_segment(aes(x = p, xend = p,
                                     y = y.start, yend = statsAdj[p]), alpha = .75, color = colors,
                                 size = 0.5) +
           scale_x_continuous(limits = c(0,length(statsAdj)), expand = c(0, 0)) +
           scale_y_continuous(limits = c(-1,1), expand = c(0, 0)) + xlab(NULL) + ylab(NULL) + 
           theme(panel.background = element_blank(), 
                 axis.text = element_blank(),
                 axis.ticks = element_blank(), 
                 panel.grid = element_blank(), 
                 axis.title = element_blank(), 
                 axis.line = element_blank(), 
                 panel.margin = element_blank(), 
                 plot.margin = rep(unit(0, "null"), 4), panel.margin = rep(unit(0, 
                                                                                "null"), 4)),
         textGrob(sprintf("%.3f", annotation$NES)), 
         textGrob(sprintf("%.2e", annotation$pval)),
         textGrob(sprintf("%.2e",annotation$padj)),
         # textGrob(DescTools::StrTrunc(paste0(names(statsAdj)[p.DEG], collapse = ' '), maxlen = 45)),
         textGrob(ifelse(nchar(paste0(names(statsAdj)[p.DEG], collapse = ' ') )>45,
                          paste0(DescTools::StrTrunc(paste0(names(statsAdj)[p.DEG], collapse = ' '), maxlen = 43), paste0('(',length(p.DEG), ')')),
                          paste0(paste0(names(statsAdj)[p.DEG], collapse = ' '))) ),
         # textGrob(DescTools::StrTrunc(paste0(names(statsAdj)[p.Immu], collapse = ' '), maxlen = 45))
         textGrob(ifelse(nchar(paste0(names(statsAdj)[p.Immu], collapse = ' ') )>45,
                          paste0(DescTools::StrTrunc(paste0(names(statsAdj)[p.Immu], collapse = ' '), maxlen = 43), paste0('(',length(p.Immu), ')')),
                          paste0(paste0(names(statsAdj)[p.Immu], collapse = ' '))) )
      )
  })
  rankPlot <- ggplot() + geom_blank() + 
    scale_x_continuous(limits = c(0,length(statsAdj)), expand = c(0, 0)) + 
    scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) + 
    xlab(NULL) + ylab(NULL) +
    theme(panel.background = element_blank(),
          axis.line = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          plot.margin = unit(c(0, 0, 0.5, 0), "npc"), 
          panel.margin = unit(c(0, 0, 0, 0), "npc"))
  
  grid.arrange(grobs = c(lapply(c("Pathways", "Ranked genes", 
                                  "NES", "pval", "padj", "DEGs", 'InnateImmune'), textGrob), unlist(ps, recursive = FALSE), 
                         list(nullGrob(), rankPlot)), 
               ncol = 7, widths = colwidths,
               top = textGrob(title, ))
}



get_run <- function (run) {
  if (run == 1) {
    print('run all')
    which_virus = c('REO', 'VSV', 'EMCV', 'MVM', 'none')
    which_treatment = c('i', 'p', 'm')
  }
  if (run == 2) {
    # treatment_plus_virus_from_treatment: cells' response to RNA viruses after polyIC treatment
    # treatment_plus_virus_from_virus: cells' response to polyIC independent of viral infection (poise).
    print('run RNA polyIC') # response expected here
    which_virus = c('REO', 'VSV', 'EMCV', 'none')
    which_treatment = c('p', 'm')
  }
  if (run == 3) {
    # treatment_plus_virus_from_treatment: cells' response to RNA viruses after interferon treatment
    # treatment_plus_virus_from_virus: cells' response to interferon independent of viral infection (poise).
    print('run RNA interfernon')
    which_virus = c('REO', 'VSV', 'EMCV', 'none')
    which_treatment = c('i', 'm')
  }
  if (run == 4) {
    print('run ssRNA interferon') # reponse expected here
    which_virus = c('VSV', 'EMCV', 'none')
    which_treatment = c('i', 'm')
  }
  if (run == 5) {
    print('run ssRNA polyIC')
    which_virus = c('VSV', 'EMCV', 'none')
    which_treatment = c('p', 'm')
  }
  if (run == 6) {
    # treatment_plus_virus_from_treatment: cells' response to RNA viruses after polyIC treatment
    # treatment_plus_virus_from_virus: cells' response to polyIC independent of viral infection (poise).
    print('run RNA all') # response expected here
    which_virus = c('REO', 'VSV', 'EMCV', 'none')
    which_treatment = c('i', 'p', 'm')
  }
  if (run == 7) {
    print('run ssRNA all')
    which_virus = c('VSV', 'EMCV', 'none')
    which_treatment = c('i', 'p', 'm')
  }
  if (run == 8) {
    print('run DNA all')
    which_virus = c('MVM', 'none')
    which_treatment = c('i', 'p', 'm')
  }
  if (run == 9) {
    print('run DNA all')
    which_virus = c('MVM', 'none')
    which_treatment = c('p', 'm')
  }
  if (run == 10) {
    print('run DNA all')
    which_virus = c('MVM', 'none')
    which_treatment = c('i', 'm')
  }
  if (run == 11) {
    print('run all virus polyIC')
    which_virus = c('REO', 'VSV', 'EMCV', 'MVM', 'none')
    which_treatment = c('p', 'm')
  }
  if (run == 12) {
    print('run all virus polyIC')
    which_virus = c('REO', 'VSV', 'EMCV', 'MVM', 'none')
    which_treatment = c('i', 'm')
  }
  return(list(which_virus = which_virus,
              which_treatment = which_treatment))
}
