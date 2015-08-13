cntOverlaps <- function(A, B, maxgap=0L, ...){
    if(maxgap>0) {
        A <- expandGR(A, maxgap)
    }
    
    numOverlaps(A, B, ...)
}