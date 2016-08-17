tileGRanges <- function(targetRegions, windowSize, step, keepPartialWindow=FALSE, ...){
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