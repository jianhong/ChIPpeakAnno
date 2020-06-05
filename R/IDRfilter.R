#' Filter peaks by IDR (irreproducible discovery rate)
#' 
#' Using IDR to assess the consistency of replicate experiments and obtain a
#' high-confidence single set of peaks
#' 
#' 
#' @param peaksA,peaksB peaklist, \link[GenomicRanges:GRanges-class]{GRanges}
#' object.
#' @param bamfileA,bamfileB file path of bam files.
#' @param maxgap,minoverlap Used in the internal call to \code{findOverlaps()}
#' to detect overlaps. See
#' \code{?\link[IRanges:findOverlaps-methods]{findOverlaps}} in the
#' \pkg{IRanges} package for a description of these arguments.
#' @param singleEnd (Default TRUE) A logical indicating if reads are single or
#' paired-end.
#' @param IDRcutoff If the IDR no less than IDRcutoff, the peak will be
#' removed.
#' @param ...  Not used.
#' @return An object GRanges
#' @author Jianhong Ou
#' @references Li, Qunhua, et al. "Measuring reproducibility of high-throughput
#' experiments." The annals of applied statistics (2011): 1752-1779.
#' @keywords misc
#' @export
#' @importFrom GenomicAlignments summarizeOverlaps Union
#' @importFrom SummarizedExperiment assay
#' @examples
#' 
#'   if(interactive()){
#'     path <- system.file("extdata", "reads", package="MMDiffBamSubset")
#'     if(file.exists(path)){
#'             bamfileA <- file.path(path, "reads", "WT_2.bam")
#'             bamfileB <- file.path(path, "reads", "Resc_2.bam")
#'             WT.AB2.Peaks <- file.path(path, "peaks", "WT_2_Macs_peaks.xls")
#'             Resc.AB2.Peaks <- file.path(path, "peaks",
#'                                        "Resc_2_Macs_peaks.xls")
#'             peaksA=toGRanges(WT.AB2.Peaks, format="MACS")
#'             peaksB=toGRanges(Resc.AB2.Peaks, format="MACS")
#'             library(idr)
#'             library(DelayedArray)
#'             IDRfilter(peaksA, peaksB,
#'                         bamfileA, bamfileB)
#'     }
#'   }
#' 
IDRfilter <- function(peaksA, peaksB, bamfileA, bamfileB, 
                      maxgap=-1L, minoverlap=0L, singleEnd=TRUE,
                      IDRcutoff=0.01, ...){
    if(!requireNamespace("idr", quietly = TRUE)){
        stop("The 'idr' package is required.")
    }
    if(!requireNamespace("DelayedArray", quietly = TRUE)){
        stop("The 'DelayedArray' package is required.")
    }
    stopifnot(is(peaksA,"GRanges"))
    stopifnot(is(peaksB,"GRanges"))
    stopifnot(file.exists(bamfileA))
    stopifnot(file.exists(bamfileB))
    ol <- findOverlapsOfPeaks(peaksA, peaksB, 
                              maxgap=maxgap, 
                              minoverlap=minoverlap)
    ol <- ol$peaklist[grepl("\\/\\/\\/", names(ol$peaklist))][[1]]
    if(length(ol)<1) return(GRanges())
    names(ol) <- paste0("olp", 
                        formatC(1:length(ol), 
                                width=nchar(as.character(length(ol))),
                                flag="0"))
    coverage <- summarizeOverlaps(features = ol, 
                                  reads = c(bamfileA, bamfileB),
                                  mode=Union, 
                                  ignore.strand = FALSE, 
                                  singleEnd=singleEnd)
    idr <- idr::est.IDR(assay(coverage)/width(
        DelayedArray::rowRanges(coverage)),
                   mu=2.07, sigma=1.34, rho=0.89, p=0.84)
    ol[idr$IDR<IDRcutoff]
}
