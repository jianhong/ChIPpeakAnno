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
    if(class(cvglists)!="list"){
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
    feature_type <- c("upstream", "5UTR", "CDS", "3UTR", "downstream")
    features <- features[features$feature_type %in% c("5UTR", "CDS", "3UTR")]
    txs <- unique(features$tx_name)
    txs <- txs[!is.na(txs)]
    if(includeIntron){
        cds <- toGRanges(TxDb, feature="CDS")
        tx_name <- cds$tx_name
        l <- lengths(tx_name)
        cds <- cds[rep(seq_along(cds), l)]
        cds$tx_name <- unlist(tx_name)
        cds$gene_id <- NULL
        cds$feature_type <- "CDS"
        txs <- intersect(cds$tx_name, txs)
        cds <- cds[cds$tx_name %in% txs]
        features <- features[features$tx_name %in% txs]
        features <- features[features$feature_type %in% c("5UTR", "3UTR")]
        features <- c(features, cds)
    }
    ## resort by txs
    features <- features[order(as.numeric(factor(features$tx_name, 
                                                 levels = txs)))]
    ## add upstream and downstream
    if(upstream.cutoff>0 && downstream.cutoff>0){
        transcripts <- toGRanges(TxDb, feature="transcript")
        transcripts <- transcripts[transcripts$tx_name %in% txs]
        suppressWarnings({
            upstream <- promoters(transcripts,
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
        transcripts.rev.strand <- revert.strand(transcripts)
        suppressWarnings({downstream <- 
            promoters(transcripts.rev.strand,
                      upstream = downstream.cutoff,
                      downstream = 0)})
        downstream <- revert.strand(downstream)
        downstream <- trim(downstream)
        downstream <- downstream[width(downstream)>=nbinsDownstream]
        upstream$gene_id <- NULL
        downstream$gene_id <- NULL
        upstream$feature_type <- "upstream"
        downstream$feature_type <- "downstream"
        features <- c(features, upstream, downstream)
        features <- 
            features[order(as.numeric(factor(features$tx_name, 
                                             levels = txs)),
                           as.numeric(factor(features$feature_type,
                                             levels = feature_type)),
                           ifelse(strand(features)=="-", -1, 1)*start(features))]
    }else{
        ## make sure features are sorted by 5' --> 3'
        features <- 
            features[order(as.numeric(factor(features$tx_name, 
                                             levels = txs)),
                           as.numeric(factor(features$feature_type,
                                             levels = c("5UTR", "CDS", "3UTR"))),
                           ifelse(strand(features)=="-", -1, 1)*start(features))]
    }
    ## split the features by seqnames and generate Views for cvglist
    seqn <- Reduce(intersect, lapply(cvglists, names))
    seqn <- intersect(seqn, seqlevels(features))
    if(length(seqn)<1){
        stop("Please check the names of cvglist. ",
             "None of them in the seqlevels of TxDb.")
    }
    features <- features[seqnames(features) %in% seqn]
    ## make sure all the region involved in same transcript should have same coverage
    features.disjoin <- disjoin(features, with.revmap=TRUE, ignore.strand=TRUE)
    revmap.len <- lengths(features.disjoin$revmap)
    revmap.id <- rep(seq_along(features.disjoin), revmap.len)
    revmap.len <- rep(revmap.len, revmap.len)
    revmap.len <- data.frame(id=unlist(features.disjoin$revmap), 
                             id2=revmap.id,
                             len=revmap.len)
    revmap.len$tx_name <- features[revmap.len$id]$tx_name
    revmap.len$feature_type <- features[revmap.len$id]$feature_type
    revmap.len$tx_type <- paste(revmap.len$tx_name, 
                                revmap.len$feature_type, sep="_")
    revmap.len$tx_type_len <- paste(revmap.len$tx_type, 
                                    revmap.len$len, sep="_")
    revmap.len$width <- width(features.disjoin[revmap.len$id2])
    revmap.width <- rowsum(revmap.len$width, revmap.len$tx_type_len)
    revmap.len$width2 <- revmap.width[revmap.len$tx_type_len, 1]
    maxWidth <- aggregate(revmap.len$width2, by=list(id=revmap.len$tx_type), 
                          FUN=max)
    revmap.len <- revmap.len[paste(revmap.len$tx_type, revmap.len$width2) %in% 
                                 paste(as.character(maxWidth[, 1]), 
                                       maxWidth[, 2]), ]
    revmap.len <- revmap.len[order(revmap.len$id), ]
    features <- features[revmap.len$id]
    ranges(features) <- ranges(features.disjoin[revmap.len$id2])
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
    ## calculate fators to balance the region be count multiple times
    len <- length(unique(features$tx_name))
    features.s <- split(features, seqnames(features))
    features.s <- features.s[seqn]
    cvglists <- lapply(cvglists, function(.ele) .ele[seqn])
    features.l <- as(features.s, "RangesList")
    cvglists.sub <- lapply(cvglists, function(.ele) .ele[features.l])
    resetGRanges <- function(ir, bins=c("upstream"=nbinsUpstream,
                                        "5UTR"=nbinsUTR,
                                        "CDS" =nbinsCDS,
                                        "3UTR"=nbinsUTR,
                                        "downstream"=nbinsDownstream)){
        stopifnot(is(ir, "IRanges"))
        if(length(ir)<1) return(ir)
        w <- c(0, width(ir))
        s <- cumsum(w)
        x <- IRanges(start = s[-length(s)]+1, end = s[-1])
        mcols(x) <- mcols(ir)
        ## merge 5UTR, CDS, and 3UTR
        mcols(x)$oid <- seq_along(ir)
        x <- split(x, mcols(x)$feature_type)
        x <- lapply(x, function(.ele){
            w <- rowsum(width(.ele), mcols(.ele)$tx_name)
            .ele <- .ele[!duplicated(mcols(.ele)$tx_name)]
            width(.ele) <- w[mcols(.ele)$tx_name, 1]
            .ele
        })
        ## split Feture By Bins
        x <- mapply(function(.x, n){
            mc <- mcols(.x)
            .x <- IRanges::tile(x=.x, n=n)
            #all(lengths(.x)==n)
            .x <- unlist(.x)
            nc <- nrow(mc)
            mc <- mc[rep(seq.int(nc), each=n), ]
            mc$id <- rep(seq.int(n), times=nc)
            mcols(.x) <- mc
            .x
        }, x, bins[names(x)], SIMPLIFY = FALSE)
        x <- unlist(RangesList(x), use.names = FALSE)
        x <- x[order(mcols(x)$oid)]
        mcols(x)$oid <- NULL
        x
    }
    features.n <- IRangesList(lapply(features.l, resetGRanges))
    cnts <- lapply(cvglists.sub, function(cvglist.sub){
        vw <- Views(cvglist.sub, features.n)
        cnt <- lapply(vw, viewMeans)
        unlist(cnt, use.names = FALSE)
    })
    features.n <- unlist(features.n)
    d <- as.data.frame(mcols(features.n)[, -1])
    d <- cbind(d, do.call(cbind, cnts))
    d <- split(d, d$feature_type)
    d <- lapply(d, function(.ele) rowsum(.ele[, -(1:2)], .ele$id, 
                                         reorder = FALSE)/len)
    d <- d[feature_type[feature_type %in% names(d)]]
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
    dat <- lapply(dat, function(.ele) rbind(.ele, NA))
    bins <- elementNROWS(dat)
    dat <- do.call(rbind, dat)
    matplot(dat, type = "l", xaxt="n", ...)
    bins.sum <- cumsum(bins)
    bins.sum <- bins.sum[-length(bins.sum)]
    abline(v = bins.sum, lty=2)
    axis(side = 1, at = bins.sum, 
         labels = rep("", length(bins.sum)), 
         tick = TRUE, lwd=-1, lwd.ticks=1)
    bins.sum2 <- c(0, bins.sum) + bins/2
    axis(side = 1, at = bins.sum2, 
         labels = feature_type, lwd=-1)
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
    legend.arg <- c(..., legend.arg)
    legend.arg <- legend.arg[!duplicated(names(legend.arg))]
    legend.arg <- legend.arg[legend.arg.names]
    legend.arg <- legend.arg[lengths(legend.arg)>0]
    do.call(legend, legend.arg)
}

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
    if(class(cvglists)!="list"){
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
                           ifelse(strand(features)=="-", -1, 1)*start(features))]
    }else{
        ## make sure features are sorted by 5' --> 3'
        gene_id <- unique(features$gene_id)
        features <- 
            features[order(as.numeric(factor(features$gene_id, 
                                             levels = gene_id)),
                           ifelse(strand(features)=="-", -1, 1)*start(features))]
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
    features.l <- as(features.s, "RangesList")
    cvglists.sub <- lapply(cvglists, function(.ele) .ele[features.l])
    resetGRanges <- function(ir, bins=c("upstream"=nbinsUpstream,
                                        "gene"=nbinsGene,
                                        "downstream"=nbinsDownstream)){
        stopifnot(is(ir, "IRanges"))
        if(length(ir)<1) return(ir)
        w <- c(0, width(ir))
        s <- cumsum(w)
        x <- IRanges(start = s[-length(s)]+1, end = s[-1])
        mcols(x) <- mcols(ir)
        mcols(x)$oid <- seq_along(ir)
        if(length(mcols(x)$feature_type)==0) mcols(x)$feature_type <- "gene"
        x <- split(x, mcols(x)$feature_type)
        x <- lapply(x, function(.ele){
            w <- rowsum(width(.ele), mcols(.ele)$gene_id)
            .ele <- .ele[!duplicated(mcols(.ele)$gene_id)]
            width(.ele) <- w[mcols(.ele)$gene_id, 1]
            .ele
        })
        ## split Feture By Bins
        x <- mapply(function(.x, n){
            mc <- mcols(.x)
            .x <- IRanges::tile(x=.x, n=n)
            #all(lengths(.x)==n)
            .x <- unlist(.x)
            nc <- nrow(mc)
            mc <- mc[rep(seq.int(nc), each=n), ]
            mc$id <- rep(seq.int(n), times=nc)
            mcols(.x) <- mc
            .x
        }, x, bins[names(x)], SIMPLIFY = FALSE)
        x <- unlist(RangesList(x), use.names = FALSE)
        x <- x[order(mcols(x)$oid)]
        mcols(x)$oid <- NULL
        x
    }
    features.n <- IRangesList(lapply(features.l, resetGRanges))
    cnts <- lapply(cvglists.sub, function(cvglist.sub){
        vw <- Views(cvglist.sub, features.n)
        cnt <- lapply(vw, viewMeans)
        unlist(cnt, use.names = FALSE)
    })
    features.n <- unlist(features.n)
    d <- as.data.frame(mcols(features.n)[, -1])
    d <- cbind(d, do.call(cbind, cnts))
    d <- split(d, d$feature_type)
    d <- lapply(d, function(.ele) rowsum(.ele[, -(1:2)], .ele$id)/len)
    feature_type <- c("upstream", "gene", "downstream")
    d <- d[feature_type[feature_type %in% names(d)]]
    return(d)
}
