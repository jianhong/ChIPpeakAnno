#' Slide windows on a given \link[GenomicRanges:GRanges-class]{GRanges} object
#' 
#' tileGRanges returns a set of genomic regions by sliding the windows in a
#' given step. Each window is called a "tile".
#' 
#' 
#' @param targetRegions A \link[GenomicRanges:GRanges-class]{GRanges} object of
#' genomic regions of interest.
#' @param windowSize Size of windows
#' @param step Step of windows
#' @param keepPartialWindow Keep last partial window or not.
#' @param ...  Not used.
#' @return A \link[GenomicRanges:GRanges-class]{GRanges} object.
#' @author Jianhong Ou
#' @keywords misc
#' @export
#' @examples
#' 
#'     genes <- GRanges(
#'         seqnames = c(rep("chr2L", 4), rep("chr2R", 5), rep("chr3L", 2)),
#'         ranges = IRanges(c(1000, 3000, 4000, 7000, 2000, 3000, 3600, 
#'                            4000, 7500, 5000, 5400), 
#'                          width=c(rep(500, 3), 600, 900, 500, 300, 900, 
#'                                  300, 500, 500),
#'                          names=letters[1:11])) 
#'     se <- tileGRanges(genes, windowSize=50, step=10)
#' 
tileGRanges <- function(targetRegions, windowSize, 
                        step, keepPartialWindow=FALSE, ...){
    if(any(width(targetRegions)<windowSize) | any(width(targetRegions)<step)){
        warning("Some of chromosomes are smaller than windowSize or step.",
                "They will be removed.")
    }
    targetRegions <- targetRegions[width(targetRegions)>=windowSize]
    targetRegions <- targetRegions[width(targetRegions)>=step]
    tileTargetRanges <- tile(x=ranges(targetRegions), width=step)
    nt <- elementNROWS(tileTargetRanges)
    tileTargetRanges.end <- rep(end(targetRegions), nt)
    tileTargetRanges <- unlist(tileTargetRanges)
    width(tileTargetRanges) <- windowSize
    tileTargetRegions <- GRanges(rep(seqnames(targetRegions), nt), 
                                 tileTargetRanges,
                                 rep(strand(targetRegions), nt),
                                 oid=rep(1:length(targetRegions), nt))
    id <- end(tileTargetRanges) > tileTargetRanges.end
    if(keepPartialWindow){
        end(tileTargetRegions[id]) <- tileTargetRanges.end[id]
    }else{
        tileTargetRegions <- 
            tileTargetRegions[!id]
    }
    return(tileTargetRegions)
}
