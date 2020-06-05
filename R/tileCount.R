#' Perform overlap queries between reads and genome by windows
#' 
#' tileCount extends
#' \link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps} by
#' providing fixed window size and step to split whole genome into windows and
#' then do queries. It will return counts in each windows.
#' 
#' 
#' @param reads A \link[GenomicRanges:GRanges-class]{GRanges},
#' \link[GenomicRanges:GRangesList-class]{GRangesList}
#' \link[GenomicAlignments:GAlignments-class]{GAlignments},
#' \link[GenomicAlignments:GAlignmentsList-class]{GAlignmentsList},
#' \link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs} or
#' \link[Rsamtools:BamFile-class]{BamFileList} object that represents the data
#' to be counted by
#' \code{\link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps}}.
#' @param genome The object from/on which to get/set the sequence information.
#' @param windowSize Size of windows
#' @param step Step of windows
#' @param keepPartialWindow Keep last partial window or not.
#' @param mode mode can be one of the pre-defined count methods. see
#' \link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps}.
#' default is countByOverlaps, alia of countOverlaps(features, reads,
#' ignore.strand=ignore.strand)
#' @param ...  Additional arguments passed to
#' \code{\link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps}}.
#' @return A
#' \link[SummarizedExperiment:RangedSummarizedExperiment-class]{RangedSummarizedExperiment}
#' object. The assays slot holds the counts, rowRanges holds the annotation
#' from genome.
#' @author Jianhong Ou
#' @keywords misc
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @importFrom GenomicAlignments summarizeOverlaps
#' @examples
#' 
#'     fls <- list.files(system.file("extdata", package="GenomicAlignments"),
#'                   recursive=TRUE, pattern="*bam$", full=TRUE)
#'     names(fls) <- basename(fls)
#'     genes <- GRanges(seqlengths = c(chr2L=7000, chr2R=10000))
#'     se <- tileCount(fls, genes, windowSize=1000, step=500)
#' 
tileCount <- function(reads, genome, windowSize=1e6, step=1e6, 
                      keepPartialWindow=FALSE,
                      mode=countByOverlaps, ...){
    targetRegions <- seqlengths(genome)
    if(any(is.na(targetRegions))) stop("Can not get seq lengths from genome.")
    targetRegions <- GRanges(names(targetRegions), 
                             IRanges(1, width=targetRegions))
    tileTargetRegions <- tileGRanges(targetRegions, windowSize, 
                                     step, keepPartialWindow)
    se <- summarizeOverlaps(features=tileTargetRegions, reads=reads, 
                            mode=mode, ...)
    return(se)
}
