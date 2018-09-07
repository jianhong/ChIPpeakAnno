#' coverage of gene body
#' 
#' @description calculate the coverage of gene body per gene 
#' per bin.
#' @param cvglists A list of \link[IRanges]{SimpleRleList} or 
#' \link[IRanges]{RleList}. It represents the coverage for samples.
#' @param TxDb An object of \code{\link[GenomicFeatures]{TxDb}}. 
#' It is used for extracting the annotations.
#' @param upstream.cutoff,downstream.cutoff cutoff length for upstream 
#' or downstream of transcript.
#' @param nbinsGene,nbinsUpstream,nbinsDownstream The number of bins
#' for gene, upstream and downstream.
#' @param includeIntron A logical value which indicates including intron or not.
#' @param minGeneLen,maxGeneLen minimal or maximal length of gene.
#' @import IRanges
#' @import GenomicRanges
#' @import GenomicFeatures
#' @import GenomeInfoDb
#' @export 
#' @author Jianhong Ou
#' @seealso \link{binOverRegions}, \link{plotBinOverRegions}
#' @examples 
#' path <- system.file("extdata", package="ChIPpeakAnno")
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(rtracklayer)
#' files <- dir(path, "bigWig")
#' if(.Platform$OS.type != "windows"){
#' cvglists <- lapply(file.path(path, files), import,
#'                    format="BigWig", as="RleList")
#' names(cvglists) <- sub(".bigWig", "", files)
#' d <- binOverGene(cvglists, TxDb.Hsapiens.UCSC.hg19.knownGene)
#' plotBinOverRegions(d)
#' }
binOverGene <- function(cvglists, TxDb, 
                        upstream.cutoff=0L, 
                        downstream.cutoff=upstream.cutoff, 
                        nbinsGene=100L,
                        nbinsUpstream=20L,
                        nbinsDownstream=nbinsUpstream,
                        includeIntron=FALSE,
                        minGeneLen=nbinsGene,
                        maxGeneLen=Inf){
  if(inherits(cvglists, c("SimpleRleList", "RleList", "CompressedRleList"))){
    cvglistsName <- substitute(deparse(cvglists))
    cvglists <- list(cvglists)
    names(cvglists) <- cvglistsName
  }
  if(!is(cvglists, "list")){
    stop("cvglists must be a list of SimpleRleList or RleList")
  }
  cls <- sapply(cvglists, class)
  if(any(!cls %in% c("SimpleRleList", "RleList", "CompressedRleList")))
    stop("cvglists must be a list of SimpleRleList or RleList")
  stopifnot(is(TxDb, "TxDb"))
  stopifnot(maxGeneLen>minGeneLen)
  stopifnot(minGeneLen>=nbinsGene)
  if(includeIntron){
    features <- toGRanges(TxDb, feature="transcript")
  }else{
    features <- toGRanges(TxDb, feature="exon")
  }
  ## reduce by gene
  GRapply <- function(X, FUN, by, ...){
    stopifnot(by!="by_id")
    stopifnot(is(X, "GRanges"))
    x <- GRanges(seqnames = mcols(X)[, by], 
                 ranges = ranges(X),
                 strand = strand(X))
    x <- FUN(x, ...)
    x$by_id <- as.character(seqnames(x))
    x <- 
      GRanges(seqnames = 
                as.character(seqnames(X[match(x$by_id, 
                                              mcols(X)[, by])])),
              ranges = ranges(x),
              strand = strand(x),
              by_id = x$by_id)
    colnames(mcols(x)) <- by
    seqinfo(x) <- seqinfo(X)[seqlevels(x)]
    x
  }
  features <- features[lengths(features$gene_id)>0]
  features$gene_id <- sapply(features$gene_id, `[`, 1)
  features <- GRapply(features, reduce, "gene_id")
  
  ## filter by minGeneLen, maxGeneLen
  features <- trim(features)
  f.width <- rowsum(width(features), features$gene_id, reorder = FALSE)
  features <- features[features$gene_id %in% 
                         rownames(f.width)[f.width[, 1]>=minGeneLen & 
                                             f.width[, 1]<maxGeneLen]]
  
  if(length(features)<2){
    stop("less than 2 gene left.")
  }
  ## add upstream and downstream
  if(upstream.cutoff>0 && downstream.cutoff>0){
    genes <- GRapply(features, range, "gene_id")
    suppressWarnings({
      upstream <- promoters(genes,
                            upstream = upstream.cutoff,
                            downstream = 0)})
    upstream <- trim(upstream)
    upstream <- upstream[width(upstream)>=nbinsUpstream]
    revert.strand <- function(gr){
      l <- levels(strand(gr))
      l.plus <- which(l=="+")
      l.minus <- which(l=="-")
      if(length(l.plus)==1){
        l[l.plus] <- "-"
      }
      if(length(l.minus)==1){
        l[l.minus] <- "+"
      }
      levels(strand(gr)) <- l
      gr
    }
    genes.rev.strand <- revert.strand(genes)
    suppressWarnings({downstream <- 
      promoters(genes.rev.strand,
                upstream = downstream.cutoff,
                downstream = 0)})
    downstream <- revert.strand(downstream)
    downstream <- trim(downstream)
    downstream <- downstream[width(downstream)>=nbinsDownstream]
    upstream$feature_type <- "upstream"
    downstream$feature_type <- "downstream"
    features$feature_type <- "gene"
    gene_id <- unique(features$gene_id)
    features <- c(features, upstream, downstream)
    features <- 
      features[order(as.numeric(factor(features$gene_id, 
                                       levels = gene_id)),
                     as.numeric(factor(features$feature_type,
                                       levels = c("upstream", 
                                                  "gene", 
                                                  "downstream"))),
                     start(features))]
  }else{
    ## make sure features are sorted by position
    gene_id <- unique(features$gene_id)
    features$feature_type <- "gene"
    features <- 
      features[order(as.numeric(factor(features$gene_id, 
                                       levels = gene_id)),
                     start(features))]
  }
  ## split the features by seqnames and generate Views for cvglist
  seqn <- Reduce(intersect, lapply(cvglists, names))
  seqn <- intersect(seqn, seqlevels(features))
  if(length(seqn)<1){
    stop("Please check the names of cvglist. ",
         "None of them in the seqlevels of TxDb.")
  }
  features <- features[seqnames(features) %in% seqn]
  len <- length(unique(features$gene_id))
  features.s <- split(features, seqnames(features))
  features.s <- features.s[seqn]
  cvglists <- lapply(cvglists, function(.ele) .ele[seqn])
  features.l <- as(features.s, "IntegerRangesList")
  
  bins=c("upstream"=nbinsUpstream,
         "gene"=nbinsGene,
         "downstream"=nbinsDownstream)
  
  features.view <- lapply(cvglists, function(.ele) {
    vw <- Views(.ele, features.l)
    cntByChr <- lapply(vw, function(.vw){
      .df <- elementMetadata(.vw)
      if(length(.df$feature_type) == 0) .df$feature_type <- "gene"
      .gr <- GRanges(paste(.df$feature_type, .df$gene_id), ranges(.vw), 
                     feature_type=.df$feature_type, gene_id=.df$gene_id)
      .gr.s <- split(.gr, seqnames(.gr))
      .cvg <- subject(.vw)
      .cvgs <- rep(list(.cvg), length(.gr.s))
      .gr.l <- as(.gr.s, "IntegerRangesList")
      names(.cvgs) <- names(.gr.l)
      .cvgs <- as(.cvgs, "SimpleRleList")
      .cvg.sub <- .cvgs[.gr.l]
      ### split the RleList into bins
      .ir <- IRanges(1, lengths(.cvg.sub))
      .ir <- IRanges::tile(.ir, n=bins[sub("^(.*?) .*$", "\\1", names(.cvg.sub))])
      names(.ir) <- names(.cvg.sub)
      .cnt <- viewMeans(Views(.cvg.sub, .ir))
      .cnt <- split(.cnt, sub("^(.*?) .*$", "\\1", names(.cnt)))
      ## make sure 5'->3'
      .cnt <- lapply(.cnt, function(x){
        strd <- as.character(strand(features))[match(sub("^(.*?) (.*$)", "\\2", names(x)), 
                                                    features$gene_id)]
        x <- split(x, strd)
        x <- lapply(x, do.call, what=rbind)
        if("-" %in% names(x)){
          x[["-"]] <- x[["-"]][, ncol(x[["-"]]):1, drop=FALSE]
        }
        do.call(rbind, x)
      })
      lapply(.cnt, colSums)
    })
    cntByFeature <- swapList(cntByChr)
    cntByFeature <- lapply(cntByFeature, function(.ele){
      colSums(do.call(rbind, .ele))/len
    })
  })
  d <- swapList(features.view)
  feature_type <- c("upstream", "gene", "downstream")
  d <- d[feature_type[feature_type %in% names(d)]]
  d <- lapply(d, do.call, what=cbind)
  return(d)
}
