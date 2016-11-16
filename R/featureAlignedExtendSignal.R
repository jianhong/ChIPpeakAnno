featureAlignedExtendSignal <- function(bamfiles, index=bamfiles, 
                                       feature.gr, 
                                       upstream, downstream, 
                                       n.tile=100, 
                                       fragmentLength,
                                       librarySize,
                                       pe=c("auto", "PE", "SE"),
                                       adjustFragmentLength,
                                       ...){
    message("The signal is being calculated for DNA-seq.")
    if(missing(fragmentLength)){
        stop("fragmentLength is missing")
    }
    if(!missing(adjustFragmentLength)){
      stopifnot(inherits(adjustFragmentLength, c("numeric", "integer")))
      stopifnot(length(adjustFragmentLength)==1)
    }
    stopifnot(class(bamfiles)=="character")
    stopifnot(class(feature.gr)=="GRanges")
    if(missing(upstream) | missing(downstream)){
        stop("upstream and downstream is missing")
    }
    stopifnot(inherits(upstream, c("numeric", "integer")))
    stopifnot(inherits(downstream, c("numeric", "integer")))
    stopifnot(inherits(n.tile, c("numeric", "integer")))
    stopifnot(inherits(fragmentLength, c("numeric", "integer")))
    stopifnot(inherits(librarySize, c("numeric", "integer")))
    pe <- match.arg(pe)
    upstream <- as.integer(upstream)
    downstream <- as.integer(downstream)
    n.tile <- as.integer(n.tile)
    fragmentLength <- as.integer(fragmentLength)
    librarySize <- as.integer(librarySize)
    stopifnot(all(width(feature.gr)==1))
    stopifnot(length(bamfiles)==length(index))
    
    totalBPinBin <- floor((upstream + downstream)/n.tile)# * length(feature.gr)
    feature.gr$oid <- 1:length(feature.gr)
    feature.gr.expand <- 
        suppressWarnings(promoters(feature.gr, upstream=upstream+max(fragmentLength), 
                  downstream=downstream+max(fragmentLength)))
    strand(feature.gr.expand) <- "*"
    feature.gr <- suppressWarnings(promoters(feature.gr,
                            upstream=upstream, 
                            downstream=downstream))
    
    grL <- tile(feature.gr, n=n.tile)
    grL <- lapply(grL, function(.gr){
        if(as.character(strand(.gr))[1]=="-"){
            .gr <- rev(.gr)
        }
        .gr$gpid <- 1:length(.gr)
        .gr
    })
    grL.eleLen <- elementNROWS(grL)
    grs <- unlist(GRangesList(grL))
    grs$oid <- rep(feature.gr$oid, grL.eleLen)
    
    if(pe=="auto") {
      pe <- mapply(function(file, id) 
        suppressMessages(testPairedEndBam(file, index=id)), 
                 bamfiles, index)
    }else{
      if(pe=="PE"){
        pe <- TRUE
      }else{
        pe <- FALSE
      }
    }
    param <- ScanBamParam(which=reduce(feature.gr.expand), 
                          flag=scanBamFlag(isSecondaryAlignment=FALSE,
                                          isNotPassingQualityControls=FALSE),
                          what=scanBamWhat())
    paramp <- ScanBamParam(which=reduce(feature.gr.expand), 
                           flag=scanBamFlag(isProperPair=TRUE,
                                            isSecondaryAlignment=FALSE,
                                            isNotPassingQualityControls=FALSE),
                           what=scanBamWhat())
    bams.gr <- mapply(function(f, i, p, .fLen) {
        if(!p){
            .ele <- granges(readGAlignments(f, index=i, param=param))
            start(.ele[strand(.ele)=="-"]) <- 
                end(.ele[strand(.ele)=="-"]) - .fLen + 1
            width(.ele) <- .fLen
            .ele
        }else{
            .ele <- readGAlignmentPairs(f, index=i, param=paramp)
            .ele <- granges(.ele)
            .ele
        }
    }, bamfiles, index, pe, fragmentLength, SIMPLIFY=FALSE)
    if(!missing(adjustFragmentLength)){
      bams.gr <- lapply(bams.gr, reCenterPeaks, width=adjustFragmentLength)
      fragmentLength <- adjustFragmentLength
    }
    ## count overlaps
    co <- lapply(bams.gr, countOverlaps, 
                 query=grs, ignore.strand=TRUE)
    
    countTable <- do.call(cbind, co)
    colnames(countTable) <- names(bams.gr)
    stopifnot(nrow(countTable)==length(grs))
#    sumByGpid <- rowsum(countTable, grs$gpid, reorder = FALSE)
#    signal <- t(t(sumByGpid)*1e8/librarySize*100/fragmentLength)/totalBPinBin
    countTable.list <- as.list(as.data.frame(countTable))
    signal <- mapply(function(.ele, .libsize, .fLen){
        tbl <- matrix(nrow=n.tile, ncol=length(feature.gr))
        tbl[(grs$oid-1)*n.tile + grs$gpid] <- .ele
        t(tbl)*1e8/.libsize*100/.fLen/totalBPinBin
    }, countTable.list, librarySize, fragmentLength, SIMPLIFY = FALSE)
    return(signal)
}