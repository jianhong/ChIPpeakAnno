countByOverlaps <- function(features, reads,  ignore.strand, inter.feature) {
    ## NOT work for parallel 
    countOverlaps(features, reads, ignore.strand=ignore.strand)
}



#' Perform overlap queries between reads and genomic features by bins
#' 
#' summarizeOverlapsByBins extends
#' \link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps} by
#' providing fixed window size and step to split each feature into bins and
#' then do queries. It will return counts by signalSummaryFUN, which applied to
#' bins in one feature, for each feature.
#' 
#' 
#' @param targetRegions A \link[GenomicRanges:GRanges-class]{GRanges} object of
#' genomic regions of interest.
#' @param reads A \link[GenomicRanges:GRanges-class]{GRanges},
#' \link[GenomicRanges:GRangesList-class]{GRangesList}
#' \link[GenomicAlignments:GAlignments-class]{GAlignments},
#' \link[GenomicAlignments:GAlignmentsList-class]{GAlignmentsList},
#' \link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs} or
#' \link[Rsamtools:BamFile-class]{BamFileList} object that represents the data
#' to be counted by
#' \code{\link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps}}.
#' @param windowSize Size of windows
#' @param step Step of windows
#' @param signalSummaryFUN function, which will be applied to the bins in each
#' feature.
#' @param mode mode can be one of the pre-defined count methods. see
#' \link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps}.
#' default is countByOverlaps, alia of countOverlaps(features, reads,
#' ignore.strand=ignore.strand)
#' @param ...  Additional arguments passed to
#' \code{\link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps}}.
#' @return A
#' \link[SummarizedExperiment:RangedSummarizedExperiment-class]{RangedSummarizedExperiment}
#' object. The assays slot holds the counts, rowRanges holds the annotation
#' from features.
#' @author Jianhong Ou
#' @keywords misc
#' @export
#' @importFrom GenomicAlignments summarizeOverlaps
#' @importFrom SummarizedExperiment colData SummarizedExperiment
#' @importFrom S4Vectors SimpleList aggregate
#' @examples
#' 
#'     fls <- list.files(system.file("extdata", package="GenomicAlignments"),
#'                   recursive=TRUE, pattern="*bam$", full=TRUE)
#'     names(fls) <- basename(fls)
#'     genes <- GRanges(
#'         seqnames = c(rep("chr2L", 4), rep("chr2R", 5), rep("chr3L", 2)),
#'         ranges = IRanges(c(1000, 3000, 4000, 7000, 2000, 3000, 3600, 
#'                            4000, 7500, 5000, 5400), 
#'                          width=c(rep(500, 3), 600, 900, 500, 300, 900, 
#'                                  300, 500, 500),
#'                          names=letters[1:11])) 
#'     se <- summarizeOverlapsByBins(genes, fls, windowSize=50, step=10)
#' 
summarizeOverlapsByBins <- function(targetRegions, reads, 
                                    windowSize=50, step=10, 
                                    signalSummaryFUN=max, 
                                    mode=countByOverlaps, ...){
    stopifnot(is(targetRegions, "GRanges"))
    stopifnot(is.function(signalSummaryFUN))
    stopifnot(length(reads)>1)
    if(length(names(targetRegions))==0 ||
       any(duplicated(names(targetRegions)))){
        stop("duplicated or null targetRegions names.")
    }
    checkFun <- signalSummaryFUN(1:10)
    if(length(checkFun)!=1){
        stop("the output of signalSummaryFUN must be a vector with length 1")
    }
    if(!inherits(checkFun, c("numeric", "integer"))){
        stop("the output of signalSummaryFUN must be a numeric.")
    }
    ## change the targetRegions by windowSize, step
#     if(any(width(targetRegions)<windowSize) | any(width(targetRegions)<step)){
#         warning("Some of targetRegions are smaller than windowSize or step.",
#                 "They will be removed.")
#     }
#     targetRegions <- targetRegions[width(targetRegions)>=windowSize]
#     targetRegions <- targetRegions[width(targetRegions)>=step]
#     tileTargetRanges <- tile(x=ranges(targetRegions), width=step)
#     nt <- elementNROWS(tileTargetRanges)
#     tileTargetRanges.end <- rep(end(targetRegions), nt)
#     tileTargetRanges <- unlist(tileTargetRanges)
#     width(tileTargetRanges) <- windowSize
#     tileTargetRegions <- GRanges(rep(seqnames(targetRegions), nt), 
#                                  tileTargetRanges,
#                                  rep(strand(targetRegions), nt),
#                                  oid=rep(1:length(targetRegions), nt))
#     tileTargetRegions <- tileTargetRegions[end(tileTargetRanges) <= 
#                                            tileTargetRanges.end]
    tileTargetRegions <- tileGRanges(targetRegions, windowSize, 
                                     step, keepPartialWindow=FALSE)
    se <- summarizeOverlaps(features=tileTargetRegions, reads=reads, 
                            mode=mode, ...)
    cnts <- aggregate(x=assay(se), 
                      by=list(oid_USED_BY_SE_OU=tileTargetRegions$oid), 
                      FUN=signalSummaryFUN, drop=FALSE)
    se.rowRanges <- targetRegions[cnts$oid_USED_BY_SE_OU]
    rownames(cnts) <- names(se.rowRanges)
    cnts$oid_USED_BY_SE_OU <- NULL
    cnts <- as.matrix(cnts)
    SummarizedExperiment(assays=SimpleList(counts=cnts), 
                         rowRanges=se.rowRanges,
                         colData=colData(se))
}
