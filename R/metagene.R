#' peak distance to features
#' 
#' Bar plot for distance to features
#' 
#' @details the bar heatmap is indicates the peaks around features. 
#' @param peaks peak list, \link[GenomicRanges:GRanges-class]{GRanges} object or
#' a \link[GenomicRanges:GRangesList-class]{GRangesList}.
#' @param AnnotationData A \link[GenomicRanges:GRanges-class]{GRanges} object 
#' or a \link[GenomicFeatures:TxDb-class]{TxDb} object.
#' @param PeakLocForDistance Specify the location of peak for calculating
#' distance,i.e., middle means using middle of the peak to calculate distance
#' to feature, start means using start of the peak to calculate the distance to
#' feature. To be compatible with previous version, by default using start
#' @param FeatureLocForDistance Specify the location of feature for calculating
#' distance,i.e., middle means using middle of the feature to calculate
#' distance of peak to feature, TSS means using start of feature when
#' feature is on plus strand and using end of feature when feature is on minus
#' strand, geneEnd means using end of feature when feature is on plus strand
#' and using start of feature when feature is on minus strand. 
#' @param upstream,downstream numeric(1). Upstream or downstream region of
#' features to plot.
#' @export
#' @importFrom ggplot2 geom_bin2d ggplot theme_bw aes_string
#' @examples 
#' path <- system.file("extdata", package="ChIPpeakAnno")
#' files <- dir(path, "broadPeak")
#' peaks <- sapply(file.path(path, files), toGRanges, format="broadPeak")
#' peaks <- GRangesList(peaks)
#' names(peaks) <- sub(".broadPeak", "", basename(names(peaks)))
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' metagenePlot(peaks, TxDb.Hsapiens.UCSC.hg19.knownGene)
metagenePlot <- function(peaks, AnnotationData, 
                     PeakLocForDistance = c("middle", "start", "end"),
                     FeatureLocForDistance = c("TSS", "middle", "geneEnd"),
                     upstream=100000, 
                     downstream=100000){
  stopifnot("peaks must be an object of GRanges or GRangesList"=
              inherits(peaks, c("GRanges", "GRangesList")))
  FeatureLocForDistance <- match.arg(FeatureLocForDistance)
  PeakLocForDistance <- match.arg(PeakLocForDistance)
  if(is(peaks, "GRanges")){
    n <- deparse(substitute(peaks))
    peaks <- GRangesList(peaks)
    names(peaks) <- n
    isGRanges <- TRUE
  }else{
    isGRanges <- FALSE
  }
  stopifnot("AnnotationData must be an object of TxDb"=
              inherits(AnnotationData, c("TxDb", "GRanges")))
  if(is(AnnotationData, "TxDb")){
    suppressMessages(g <- genes(AnnotationData))
  }else{
    g <- AnnotationData
    if(length(names(g))!=length(g)){
      names(g) <- paste0("ann", 
                         formatC(seq_along(g), 
                                 width = nchar(as.character(length(g))),
                                 flag = 0))
    }
  }
  seql <- seqlevelsStyle(g)
  seqn <- seqlevels(g)
  peaks <- lapply(peaks, function(.ele){
    seqlevelsStyle(.ele) <- seql[1]
    .ele <- .ele[seqnames(.ele) %in% seqn]
    seqlevels(.ele) <- seqn
    .ele
  })
  features <- switch(FeatureLocForDistance,
                     "TSS"=promoters(g, upstream=0, downstream=1),
                     "geneEnd"=downstreams(g, upstream=1, downstream=0),
                     "middle"=reCenterPeaks(g, width=1))
  suppressWarnings({
    features.ext <- switch(FeatureLocForDistance,
                         "TSS"=promoters(g, upstream=upstream, 
                                         downstream=downstream),
                         "geneEnd"=downstreams(g, upstream=upstream, 
                                           downstream=downstream),
                         "middle"=reCenterPeaks(g, width=upstream+downstream))})
  features.ext <- GenomicRanges::trim(features.ext)
  ol <- lapply(peaks, function(.peaks){
    findOverlaps(.peaks, features.ext)
  })
  pp <- mapply(ol, peaks, FUN=function(.ol, .peaks){
    a <- .peaks[queryHits(.ol)]
    a <- switch(PeakLocForDistance,
                "start"=promoters(a, upstream=0, downstream=1),
                "end"=downstreams(a, upstream=1, downstream=0),
                "middle"=reCenterPeaks(a, width=1))
    b <- features[subjectHits(.ol)]
    d <- distance(a, b)
    s <- pcompare(ranges(a), ranges(b))
    s <- ifelse(as.character(strand(b)) == "-", -1*s, s)
    d <- data.frame(id=names(b), distance=d*sign(s))
  }, SIMPLIFY = FALSE)
  dat <- do.call(rbind, pp)
  dat$peaks <- rep(names(peaks), vapply(pp, nrow, FUN.VALUE = 0))
  ggplot(dat, aes_string(x="distance", y="peaks")) + 
    geom_bin2d()+theme_bw()
}
