countByOverlaps <- function(features, reads,  ignore.strand, inter.feature) {
    countOverlaps(features, reads, ignore.strand=ignore.strand)
}

summarizeOverlapsByBins <- function(targetRegions, reads, 
                                    windowSize=50, step=10, 
                                    signalSummaryFUN=max, 
                                    mode=countByOverlaps, ...){
    stopifnot(class(targetRegions)=="GRanges")
    stopifnot(base::mode(signalSummaryFUN)=="function")
    stopifnot(length(reads)>1)
    if(length(names(targetRegions))==0 || any(duplicated(names(targetRegions)))){
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
    if(any(width(targetRegions)<windowSize) | any(width(targetRegions)<step)){
        warning("Some of targetRegions are smaller than windowSize or step.",
                "They will be removed.")
    }
    targetRegions <- targetRegions[width(targetRegions)>=windowSize]
    targetRegions <- targetRegions[width(targetRegions)>=step]
    tileTargetRanges <- tile(x=ranges(targetRegions), width=step)
    tileTargetRanges <- mapply(function(.ele, .end){
        width(.ele) <- windowSize
        .ele[end(.ele)<=.end] ## drop the tiles of over the end
    }, tileTargetRanges, end(targetRegions))
    tileTargetRanges <- IRangesList(tileTargetRanges)
    nt <- elementNROWS(tileTargetRanges)
    tileTargetRegions <- GRanges(rep(seqnames(targetRegions), nt), 
                                 unlist(tileTargetRanges),
                                 rep(strand(targetRegions), nt),
                                 oid=rep(1:length(targetRegions), nt))
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