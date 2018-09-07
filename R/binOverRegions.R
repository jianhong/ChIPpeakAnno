#' coverage of chromosome regions
#' 
#' @description calculate the coverage of 5'UTR, CDS and 3'UTR per transcript 
#' per bin.
#' @param cvglists A list of \link[IRanges]{SimpleRleList} or 
#' \link[IRanges]{RleList}. It represents the coverage for samples.
#' @param TxDb An object of \code{\link[GenomicFeatures]{TxDb}}. 
#' It is used for extracting the annotations.
#' @param upstream.cutoff,downstream.cutoff cutoff length for upstream 
#' or downstream of transcript.
#' @param nbinsCDS,nbinsUTR,nbinsUpstream,nbinsDownstream The number of bins
#' for CDS, UTR, upstream and downstream.
#' @param includeIntron A logical value which indicates including intron or not.
#' @param minCDSLen,minUTRLen minimal length of CDS or UTR of transcript.
#' @param maxCDSLen,maxUTRLen maximal length of CDS or UTR of transctipt.
#' @import IRanges
#' @import GenomicRanges
#' @import GenomicFeatures
#' @import GenomeInfoDb
#' @export 
#' @author Jianhong Ou
#' @seealso \link{binOverGene}, \link{plotBinOverRegions}
#' @examples 
#' path <- system.file("extdata", package="ChIPpeakAnno")
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(rtracklayer)
#' files <- dir(path, "bigWig")
#' if(.Platform$OS.type != "windows"){
#' cvglists <- lapply(file.path(path, files), import,
#'                    format="BigWig", as="RleList")
#' names(cvglists) <- sub(".bigWig", "", files)
#' d <- binOverRegions(cvglists, TxDb.Hsapiens.UCSC.hg19.knownGene)
#' plotBinOverRegions(d)
#' }
binOverRegions <- function(cvglists, TxDb, 
                           upstream.cutoff=1000L, 
                           downstream.cutoff=upstream.cutoff, 
                           nbinsCDS=100L, nbinsUTR=20L, 
                           nbinsUpstream=20L,
                           nbinsDownstream=nbinsUpstream,
                           includeIntron=FALSE,
                           minCDSLen=nbinsCDS,
                           minUTRLen=nbinsUTR,
                           maxCDSLen=Inf,
                           maxUTRLen=Inf){
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
    stopifnot(maxCDSLen>minCDSLen)
    stopifnot(minCDSLen>=nbinsCDS)
    stopifnot(maxUTRLen>minUTRLen)
    stopifnot(minUTRLen>=nbinsUTR)
    features <- toGRanges(TxDb, feature="geneModel")
    features.reduced <- reduce(features)
    feature_type <- c("upstream", "5UTR", "CDS", "3UTR", "downstream")
    features <- features[features$feature_type %in% c("5UTR", "CDS", "3UTR")]
    transcripts <- transcripts(TxDb, columns=c("tx_name", "gene_id"))
    trx <- unique(mcols(transcripts))
    features$gene_id <- trx[match(features$tx_name, trx$tx_name), "gene_id"]
    features <- features[lengths(features$gene_id)>0]
    features$tx_name <- sapply(features$gene_id, function(.ele) .ele[1])
    features$gene_id <- NULL
    
    features.disjoin <- disjoin(features, with.revmap=TRUE)
    features.disjoin.1 <- features.disjoin[rep(seq_along(features.disjoin), lengths(features.disjoin$revmap))]
    features.disjoin.1$revmap <- unlist(features.disjoin$revmap)
    features.disjoin.1$feature_type <- features$feature_type[features.disjoin.1$revmap]
    features.disjoin.1$tx_name <- features$tx_name[features.disjoin.1$revmap]
    features.disjoin.1$revmap2 <- rep(seq_along(features.disjoin), lengths(features.disjoin$revmap))
    feature_type2 <- split(features.disjoin.1$feature_type, features.disjoin.1$revmap2)
    feature_type2 <- lapply(feature_type2, unique)
    feature_type2 <- feature_type2[lengths(feature_type2)==1]
    features.disjoin.1 <- features.disjoin.1[match(as.numeric(names(feature_type2)), features.disjoin.1$revmap2)]
    features.disjoin.1$seqn <- paste(features.disjoin.1$tx_name, 
                                     features.disjoin.1$feature_type, 
                                     as.character(seqnames(features.disjoin.1)))
    features1 <- GRanges(seqnames = features.disjoin.1$seqn, ranges = ranges(features.disjoin.1),
                         strand = strand(features.disjoin.1))
    features1 <- reduce(features1)
    mcols(features1) <- data.frame(do.call(rbind, strsplit(as.character(seqnames(features1)), " ")), 
                                   stringsAsFactors=FALSE)
    colnames(mcols(features1)) <- c("tx_name", "feature_type", "seqn")
    features1 <- GRanges(seqnames = features1$seqn, 
                         ranges = ranges(features1),
                         strand = strand(features1),
                         tx_name = features1$tx_name,
                         feature_type = features1$feature_type)
    seqinfo(features1) <- seqinfo(features)[seqlevels(features1)]
    features <- features1
    rm(list=c("features.disjoin", "features.disjoin.1", "feature_type2", "features1"))
    
    txs <- unique(features$tx_name)
    txs <- txs[!is.na(txs)]
    if(includeIntron){
        cds <- features[features$feature_type=="CDS"]
        cds <- GRanges(seqnames = paste(as.character(seqnames(cds)), cds$tx_name),
                       ranges = ranges(cds), strand = strand(cds))
        cds <- reduce(cds, min.gapwidth=1e9)
        mcols(cds) <- data.frame(do.call(rbind, strsplit(as.character(seqnames(cds)), " ")), 
                                 stringsAsFactors = FALSE)
        colnames(mcols(cds)) <- c("seqn", "tx_name")
        cds <- GRanges(seqnames = cds$seqn, ranges = ranges(cds), strand = strand(cds),
                       tx_name = cds$tx_name, feature_type = "CDS")
        seqinfo(cds) <- seqinfo(features)[seqlevels(cds)]
        features <- features[features$feature_type %in% c("5UTR", "3UTR")]
        features <- c(features, cds)
    }
    ## resort by txs
    features <- features[order(as.numeric(factor(features$tx_name, 
                                                 levels = txs)))]
    ## features must contain 5UTR, CDS and 3UTR
    txs <- split(features$feature_type, features$tx_name)
    txs <- lapply(txs, unique)
    txs <- names(txs)[lengths(txs)==3]
    features <- features[features$tx_name %in% txs]
    ## add upstream and downstream
    if(upstream.cutoff>0 && downstream.cutoff>0){
        genes <- GRanges(seqnames = paste(as.character(seqnames(features)), features$tx_name),
                         ranges = ranges(features), strand = strand(features))
        genes <- reduce(genes, min.gapwidth=1e9)
        mcols(genes) <- data.frame(do.call(rbind, strsplit(as.character(seqnames(genes)), " ")),
                                   stringsAsFactors = FALSE)
        colnames(mcols(genes)) <- c("seqn", "tx_name")
        genes <- GRanges(seqnames = genes$seqn, ranges = ranges(genes), strand = strand(genes),
                       tx_name = genes$tx_name)
        seqinfo(genes) <- seqinfo(features)[seqlevels(genes)]
        suppressWarnings({
            upstream <- promoters(genes,
                                  upstream = upstream.cutoff,
                                  downstream = 0)})
        upstream <- trim(upstream)
        upstream <- upstream[width(upstream)>=nbinsUpstream]
        revert.strand <- function(gr){
            l <- as.character(strand(gr))
            l.plus <- which(l=="+")
            l.minus <- which(l=="-")
            if(length(l.plus)>=1){
                l[l.plus] <- "-"
            }
            if(length(l.minus)>=1){
                l[l.minus] <- "+"
            }
            strand(gr) <- l
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
        #remove the overlapping parts
        rmOverlaps <- function(gr, tobeRemoved=features.reduced){
          gr.intersect <- intersect(gr, tobeRemoved, ignore.strand=TRUE)
          gr.intersect$tx_name <- NA
          gr.intersect$feature_type <- NA
          gr.disjoin <- disjoin(c(gr, gr.intersect), 
                                      with.revmap=TRUE, ignore.strand=TRUE)
          gr.disjoin <- gr.disjoin[lengths(gr.disjoin$revmap)==1]
          gr.disjoin$revmap <- unlist(gr.disjoin$revmap)
          strand(gr.disjoin) <- strand(gr[gr.disjoin$revmap])
          mcols(gr.disjoin) <- mcols(gr[gr.disjoin$revmap])
          gr.disjoin
        }
        upstream <- rmOverlaps(upstream)
        downstream <- rmOverlaps(downstream)
        features <- c(features, upstream, downstream)
    }
    ## split the features by seqnames and generate Views for cvglist
    seqn <- Reduce(intersect, lapply(cvglists, names))
    seqn <- intersect(seqn, seqlevels(features))
    if(length(seqn)<1){
        stop("Please check the names of cvglist. ",
             "None of them in the seqlevels of TxDb.")
    }
    features <- features[seqnames(features) %in% seqn]
    
    ## make sure no overlaps in all the region 
    features.disjoin <- disjoin(features, with.revmap=TRUE)
    features.disjoin <- features.disjoin[lengths(features.disjoin$revmap)==1]
    features.disjoin$revmap <- unlist(features.disjoin$revmap)
    features.disjoin$feature_type <- features$feature_type[features.disjoin$revmap]
    features.disjoin$tx_name <- features$tx_name[features.disjoin$revmap]
    features.disjoin$revmap <- NULL
    features.disjoin.tx_name <- split(features.disjoin$feature_type, features.disjoin$tx_name)
    features.disjoin.tx_name <- lapply(features.disjoin.tx_name, unique)
    features.disjoin.tx_name <- 
      names(features.disjoin.tx_name)[lengths(features.disjoin.tx_name) ==
                                        max(lengths(features.disjoin.tx_name))]
    features.disjoin <- features.disjoin[features.disjoin$tx_name %in% 
                                           features.disjoin.tx_name]
    seqinfo(features.disjoin) <- seqinfo(features)[seqlevels(features.disjoin)]
    features <- features.disjoin
    rm(features.disjoin)
    ## make sure features are sorted by pos
    features <- 
      features[order(as.numeric(factor(features$tx_name, 
                                       levels = txs)),
                     as.numeric(factor(features$feature_type,
                                       levels = feature_type[feature_type %in% unique(features$feature_type)])),
                     start(features))]
    
    ## filter by minCDSLen, min5UTR, min3UTR
    filterByLen <- function(type, minlen, maxLen){
        x <- features[features$feature_type %in% type]
        if(length(x)<1) return(x)
        x.width <- rowsum(width(x), x$tx_name, reorder = FALSE)
        rownames(x.width)[x.width[, 1] >= minlen & x.width[, 1] < maxLen]
    }
    txs <- Reduce(intersect, mapply(filterByLen, 
                                    c("CDS", "5UTR", "3UTR", 
                                      "upstream", "downstream"), 
                                    c(minCDSLen, minUTRLen, minUTRLen, 
                                      nbinsUpstream, nbinsDownstream),
                                    c(maxCDSLen, maxUTRLen, maxUTRLen,
                                      Inf, Inf)))
    if(length(txs)<2){
        stop("less than 2 CDS left.")
    }
    features <- features[features$tx_name %in% txs]
    ## features must contain 5UTR, CDS, 3UTR, and/or upstream, downstream.
    txs <- split(features$feature_type, features$tx_name)
    txs <- lapply(txs, unique)
    txs <- names(txs)[lengths(txs)==length(unique(features$feature_type))]
    if(length(txs)<2){
      stop("less than 2 CDS left.")
    }
    features <- features[features$tx_name %in% txs]
    ## calculate fators to balance the region be count multiple times
    len <- length(unique(features$tx_name))
    features.s <- split(features, seqnames(features))
    features.s <- features.s[seqn]
    cvglists <- lapply(cvglists, function(.ele) .ele[seqn])
    features.l <- as(features.s, "IntegerRangesList")
    
    bins=c("upstream"=nbinsUpstream,
           "5UTR"=nbinsUTR,
           "CDS" =nbinsCDS,
           "3UTR"=nbinsUTR,
           "downstream"=nbinsDownstream)
    
    features.view <- lapply(cvglists, function(.ele) {
      vw <- Views(.ele, features.l)
      cntByChr <- lapply(vw, function(.vw){
        .df <- elementMetadata(.vw)
        .gr <- GRanges(paste(.df$feature_type, .df$tx_name), ranges(.vw), 
                       feature_type=.df$feature_type, tx_name=.df$tx_name)
        .gr.s <- split(.gr, seqnames(.gr))
        .cvg <- subject(.vw)
        .cvgs <- rep(list(.cvg), length(.gr.s))
        .gr.l <- as(.gr.s, "IntegerRangesList")
        names(.cvgs) <- names(.gr.l)
        .cvgs <- as(.cvgs, "SimpleRleList")
        .cvg.sub <- .cvgs[.gr.l]
        ### split the RleList into bins
        .ir <- IRanges(rep(1, length(.cvg.sub)), lengths(.cvg.sub))
        .ir <- IRanges::tile(.ir, n=bins[sub("^(.*?) .*$", "\\1", names(.cvg.sub))])
        names(.ir) <- names(.cvg.sub)
        .cnt <- viewMeans(Views(.cvg.sub, .ir))
        .cnt <- split(.cnt, sub("^(.*?) .*$", "\\1", names(.cnt)))
        ## make sure 5'->3'
        .cnt <- lapply(.cnt, function(x){
          strd <- as.character(strand(features))[match(sub("^(.*?) (.*$)", "\\2", names(x)), 
                                                       features$tx_name)]
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
    d <- d[feature_type[feature_type %in% names(d)]]
    d <- lapply(d, do.call, what=cbind)
    return(d)
}


#' plot the coverage of regions
#' 
#' @description plot the output of \link{binOverRegiions} or \link{binOverGene}
#' @param dat A list of matrix which indicate the coverage of regions per bin
#' @param ... Parameters could be used by \link[graphics]{matplot}
#' @export
#' @author Jianhong Ou
#' @seealso \link{binOverRegions}, \link{binOverGene}
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
plotBinOverRegions <- function(dat, ...){
    feature_type <- c("upstream", "5UTR", "CDS", "gene", "3UTR", "downstream")
    if(!all(names(dat) %in% feature_type)){
        stop("names of dat must be upstream, 5UTR, CDS, 3UTR or downstream")
    }
    feature_type <- feature_type[feature_type %in% names(dat)]
    dat <- dat[feature_type]
    bins <- elementNROWS(dat)
    dat <- do.call(rbind, dat)
    matplot(dat, type = "l", xaxt="n", ...)
    bins.sum <- cumsum(bins) + .5
    bins.sum <- bins.sum[-length(bins.sum)]
    abline(v = bins.sum, lty=2, col="gray")
    axis(side = 1, at = bins.sum, 
         labels = rep("", length(bins.sum)), 
         tick = TRUE, lwd=-1, lwd.ticks=1, ...)
    bins.sum2 <- c(0, bins.sum) + bins/2
    axis(side = 1, at = bins.sum2, 
         labels = feature_type, lwd=-1, ...)
    legend.arg <- formals("legend")
    legend.arg.names <- names(legend.arg)
    legend.arg$x <- "topright"
    legend.arg$box.col <- NA
    legend.arg$merge <- FALSE
    legend.arg$legend <- colnames(dat)
    matplot.arg <- formals("matplot")
    for(arg in c("col", "lty", "lwd", "bg")){
        legend.arg[[arg]] <- matplot.arg[[arg]]
    }
    legend.arg <- c(list(...), legend.arg)
    legend.arg <- legend.arg[!duplicated(names(legend.arg))]
    legend.arg <- legend.arg[legend.arg.names]
    legend.arg <- legend.arg[lengths(legend.arg)>0]
    do.call(legend, legend.arg)
}

