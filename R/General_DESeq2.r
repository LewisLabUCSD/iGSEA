library(DESeq2)
library(gplots)

##########################
#### functions
clean_res <- function(res,ddsR){
	res$status = 'OK'
	if (any(mcols(ddsR)$cooksOutlier,na.rm=T)) res[which(mcols(ddsR)$cooksOutlier),]$status = 'Outlier'
	fil = res$baseMean > 10
	res[fil,]$padj = p.adjust(res[fil,]$pvalue, 'BH')
	res[!fil,]$padj = NA
	res[!fil,]$status = 'Low'
}

## DE function
run_deseq2<-function(expss,annot,full,reduced){
	### initialize dds object
	dds = DESeqDataSetFromMatrix(countData=expss, colData=annot, design= full)
	### run LRT
	ddsR = DESeq(dds,test='LRT',full= full , reduced = reduced)
	### get results
	results(ddsR) 
}

########################
######## code

### load data

## Annotation
annot <- read.csv(...) # annotation file ( samples x features )

## Expression
res_ex <- list()
expss <- read.csv(...,row.names=1,header=TRUE) # gene counts file ( genes x samples), gene names are 1st row

## Initialize batch effects
reduced_base = ~ batch # if there are known (e.g. plate) or computed (e.g. combat) batch effects
reduced_base = ~ 1 # if there are no batch effects considered in this analysis

## Univariate Analysis 
# feature 1 effect
reduced = reduced_base ; full = update(reduced,~ . + feature1)  ;
dds = DESeqDataSetFromMatrix(countData=expss, colData=annot, design= full)
ddsR = DESeq(dds,test='LRT',full= full , reduced = reduced)
res_ex[['feature1']] =  results(ddsR) 

# feature 2 effect
reduced = reduced_base ; full = update(reduced,~ . + feature2)  ;
dds = DESeqDataSetFromMatrix(countData=expss, colData=annot, design= full)
ddsR = DESeq(dds,test='LRT',full= full , reduced = reduced)
res_ex[['feature2']] =  results(ddsR) 

## Bivariate Analysis
# additive effect when accounting for contributions of feature1
reduced = update(reduced_base,~ . + feature1) ; full = update(reduced,~ . +feature2)  ;
dds = DESeqDataSetFromMatrix(countData=expss, colData=annot, design= full)
ddsR = DESeq(dds,test='LRT',full= full , reduced = reduced)
res_ex[['additive_wo_feature1']] =  results(ddsR) 

# additive effect when accounting for contributions of feature2
reduced = update(reduced_base,~ . + feature2) ; full = update(reduced,~ . +feature1)  ;
dds = DESeqDataSetFromMatrix(countData=expss, colData=annot, design= full)
ddsR = DESeq(dds,test='LRT',full= full , reduced = reduced)
res_ex[['additive_wo_feature2']] =  results(ddsR) 

## Interaction effect
# interaction effect when accounting for the additive effect
reduced = update(reduced_base,~ . feature1 + feature2) ; full = update(reduced,~ . + feature1:feature2)  ;
dds = DESeqDataSetFromMatrix(countData=expss, colData=annot, design= full)
ddsR = DESeq(dds,test='LRT',full= full , reduced = reduced)
res_ex[['interaction']] =  results(ddsR) 



save(res_ex,file='res_ex.rda')

#output tables
tmp=sapply(1:length(res_ex),function(i) write.csv(res_ex[[i]],file=paste('DE_out/',names(res_ex)[i],sep='')))
#output MAplots
pdf('MAplots.pdf')
tmp=sapply(1:length(res_ex),function(i) plotMA(res_ex[[i]],main=paste('MAplot:',names(res_ex[i]))))
dev.off()
# output heatmap from DE genes
# ...


# p-values
pdf('significance.pdf',width=10)
par(mfrow=c(1,2))
i=1
while(i<=length(res_ex)){
	hist(res_ex[[i]]$pvalue, main=names(res_ex)[i],xlab='Pr( chi^2(LRT) )' )
	hist(res_ex[[i]]$padj, main=names(res_ex)[i],xlab='FDR( Pr( chi^2(LRT) ))')
	i<-i+1
}
dev.off()

