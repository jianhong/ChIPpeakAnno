#' count overlaps
#' 
#' Count overlaps with max gap.
#' 
#' @param A,B A GRanges object.
#' @param maxgap A single integer >= 0.
#' @param ... parameters passed to \link[regioneR]{numOverlaps}#' 
#' @importFrom regioneR numOverlaps
cntOverlaps <- function(A, B, maxgap=0L, ...){
    if(maxgap>0) {
        A <- expandGR(A, maxgap)
    }
    
    numOverlaps(A, B, ...)
}