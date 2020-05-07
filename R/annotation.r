library(RamiGO)
library(biomaRt)
library(igraph)

#### GO hierarchy query and visualization

#' Ramigo_wrapper
#'
#' A function for querying the GO hierarchy given a vector of go terms. This function will by default (read=TRUE) return the go hierarchy as an igraph object
#' @param go a vector of GO terms to search
#' @param attempts a numeric variable indicating the number of times the server will be contacted
#' @param file a filename for the .dot file downloaded from Ramigo
#' @param picType either 'dot' or 'png' to indicate the type of file to download. 'dot' is required to read the file as a network.
#' @param read boolean variable indicating the return type: TRUE returns an igraph network descriping the hierarchy, FALSE returns the file dowloaded from Ramigo 
Ramigo_wrapper <- function(go,attempts=10,file='ramigo.dot',picType=c('dot','png'),read=TRUE,...){
	dotRes <- NULL
	while(is.null(dotRes)){
		try( dotRes <- getAmigoTree(goIDs=go, picType=picType, saveResult=TRUE,filename=file,...) )
		attempts = attempts - 1
		if(attempts==0){stop('server cannot be reached with multiple attempts')}
	}
	#tt <- readAmigoDot(object=dotRes)
	if(read){
		return(readAmigoDot_man(file)) 
	}else{
		return(dotRes)
	}

}

# internal
readAmigoDot_man <- function(file){
	r <- readLines(file)
	edges <- list()
	nodes <- list()
	for(r_i in r){
		if(grepl('->',r_i)){
			# edge
			edges[[r_i]] <- gsub('\t','',strsplit(r_i,' ')[[1]][c(1,3)])
		}else if(grepl('GO|UBERON|CARO|CL|OBA',r_i )) {
			# node
			nodes[[r_i]] <- strsplit(r_i,'>|<')[[1]][8]
		}else{
			print(paste('skip',r_i))
		}
	}
	nodes<-unlist(nodes)
	edges<-lapply(edges,function(x) unlist(nodes[as.numeric(unlist(gsub('node','',x)))]))
	eL <-do.call(rbind,edges)
	rownames(eL) <- NULL
	colnames(eL) <- c('parent','child')
	graph.data.frame(eL,directed=FALSE)
}

#### Query biomart 
# for further instruction see http://uswest.ensembl.org/info/data/biomart/biomart_r_package.html
#' biomaRt_query
#'
#' A function for query biomart for GO terms, KEGG  terms, evolutionary statistics and more!
#' @param ids vector of ids to query on biomart
#' @param id_type character choose the id_type of the variable "ids." For example filter (id_type) initialize the biomart object then runlistFilters see: ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl"); head(listFilters(ensembl))
#' @param extra vector of other datatypes to query (don't query too many, biomart will cancel the job). For example attributes (extra) initialize the biomart object then runlistFilters see: ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl"); head(listAttributes(ensembl))
#' @param biomart a string describing the archive to query: To view possible archives: listEnsembl(); listEnsembl("GRCh=37"); listEnsembl(version=78)
#' @param dataset character indicting the dataset to query: <organism>_<level>_<archive/biomart>. 
#' @param ... the archive version (version=78) or genome number can be specified (GRCh=37)
biomaRt_query <- function(ids,id_type='ensembl_gene_id',extra=NULL,biomart="ensembl",dataset="hsapiens_gene_ensembl",...){  # ENSEMBL_MART_ENSEMBL
	ensembl=useMart(biomart=biomart,dataset=dataset, host="www.ensembl.org",...)
	BM = getBM(attributes=unique(c('go_id',id_type,extra)), filters=id_type, values=ids, mart=ensembl)
	BM[BM==''] = NA  # remove missing values
	BM = na.omit(BM)
	BM
}

#### Load gtf
read.gtf<-function(gtf_file='~/scratch_ln/genomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf'){
	gtf <- read.csv(gtf_file,sep='\t',header=FALSE) #,'exons')
	tmp=strsplit(as.character(gtf$V9),';')
	gtf$gene_name = lapply(tmp,function(x) gsub('gene_name ','',x[grepl('gene_name',x)]))
	gtf$gene_id = lapply(tmp,function(x) gsub('gene_id ','',x[grepl('gene_id',x)]))
	#gtf$transcript_id = lapply(tmp,function(x) gsub('transcript_id ','',x[grepl('transcript_id',x)]))
	gtf = gtf[gtf$V3=='exon',]
	save(gtf,file='~/scratch_ln/ref_tmp.rda')

	gr_split = lapply(unique(as.character(gtf$gene_id)),function(x){
		print(x)
		gtf_i = gtf[gtf$gene_id==x,]
		reduce( tmp<-GRanges(seqnames=Rle(paste('chr',gtf_i$V1,sep='')),
			ranges=IRanges(gtf_i$V4 ,end=gtf_i$V5 ,names=gtf_i$gene_name), 
			strand=Rle(strand(gtf_i$V7)),
			type=gtf_i$V2,
			gene_name=gtf_i$gene_name,
			gene_id=gtf_i$gene_id,
			annotation=gtf_i$V9
			) )
		})
	names(gr_split) = unique(as.character(gtf$gene_id))
	len = lapply( lapply(gr_split,width ) , sum)
	save(gtf,gr_split,len,file='~/scratch_ln/ref_tmp.rda')
}