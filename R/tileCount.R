tileCount <- function(reads, genome, windowSize=1e6, step=1e6, 
                      keepPartialWindow=FALSE,
                      mode=countByOverlaps, ...){
    targetRegions <- seqlengths(genome)
    if(any(is.na(targetRegions))) stop("Can not get seq lengths from genome.")
    targetRegions <- GRanges(names(targetRegions), IRanges(1, width=targetRegions))
    tileTargetRegions <- tileGRanges(targetRegions, windowSize, 
                                     step, keepPartialWindow)
    se <- summarizeOverlaps(features=tileTargetRegions, reads=reads, 
                            mode=mode, ...)
    return(se)
}