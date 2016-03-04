##require library graph, RBGL
findOverlapsOfPeaks <- function(..., maxgap=0L, minoverlap=1L, 
                                ignore.strand=TRUE, connectedPeaks=c("min", "merge", "keepAll")){
    ###check inputs
    NAME_conn_string <- "___conn___"
    NAME_short_string <- "__"
    NAME_long_string <- "///"
    PeaksList <- list(...)
    n <- length(PeaksList)
    if(n==1){
        PeaksList <- PeaksList[[1]]
        n <- length(PeaksList)
        names <- names(PeaksList)
        if(is.null(names)) names <- paste("peaks", 1:n, sep="")
    }else{
        ##save dots arguments names
        dots <- substitute(list(...))[-1]
        names <- make.names(unlist(sapply(dots, deparse)))
    }
    if(any(grepl(NAME_short_string, names))){
        stop(paste("The name of peaks could not contain", NAME_short_string))
    }
    if(n<2){
        stop("Missing required argument Peaks!")
    }
    if(n>5){
        stop("The length of input peaks list should no more than 5")
    }
    connectedPeaks <- match.arg(connectedPeaks)
    if(any(duplicated(names)))
        stop("Same input Peaks detected!")
    
    venn_cnt <- 
        vennCounts(PeaksList, n=n, names=names, 
                   maxgap=maxgap, minoverlap=minoverlap, by="region",
                   ignore.strand=ignore.strand, connectedPeaks=connectedPeaks)
    
    outcomes <- venn_cnt$venn_cnt[, 1:n]
    xlist <- do.call(rbind, venn_cnt$xlist)
    xlist <- xlist - 1
    xlist <- xlist[nrow(xlist):1,,drop=FALSE] 
    ## reverse xlist to match the order of names
    xlist <- apply(xlist, 2, paste, collapse="")
    if(length(venn_cnt$all)!=length(xlist)) 
        stop("length of 'xlist' and 'all' should be identical.")
    all <- do.call(rbind,
                   mapply(function(.ele, .id, .gp) cbind(.id, .ele, .gp), 
                          venn_cnt$all, 1:length(venn_cnt$all), 
                          xlist, SIMPLIFY=FALSE))
    
    all.peaks <- venn_cnt$Peaks[all[,2]]
    names(all.peaks) <- gsub(NAME_conn_string,
                             NAME_short_string,
                             names(all.peaks))
    if(!is.null(all.peaks$old_strand_HH)){
        strand(all.peaks) <- all.peaks$old_strand_HH
        all.peaks$old_strand_HH <- NULL
    }
    all.peaks$gpForFindOverlapsOfPeaks <- all[, 1]
    all.peaks$gpType <- all[, 3]
    all.peaks.split <- split(all.peaks, all.peaks$gpType)
    listname <- apply(outcomes, 1, 
                      function(id) 
                          paste(names[as.logical(id)], 
                                collapse=NAME_long_string))
    listcode <- apply(outcomes, 1, paste, collapse="")
    listname <- listname[-1]
    listcode <- listcode[-1]
    names(listname) <- listcode
    names(all.peaks.split) <- listname[names(all.peaks.split)]
    peaklist <- lapply(all.peaks.split, reduce, 
                       min.gapwidth=maxgap, with.revmap=TRUE, 
                       ignore.strand=ignore.strand)
    peaklist <- mapply(function(peaks, info){
        revmap <- peaks$revmap
        peakNames <- do.call(rbind, 
                             mapply(function(id, gp) cbind(id, gp), 
                                    revmap, 1:length(revmap), 
                                    SIMPLIFY=FALSE))
        peaks$peakNames <- CharacterList(split(names(info)[peakNames[, 1]], 
                                               peakNames[, 2]), 
                                         compress=TRUE)
        if(ignore.strand){
            strand <- split(as.character(strand(info))[peakNames[, 1]], 
                            peakNames[, 2])
            strand <- lapply(strand, unique)
            l <- sapply(strand, length)
            strand[l>=2] <- "*"
            strand(peaks) <- unlist(strand)
        }
        peaks$revmap <- NULL
        peaks
    }, peaklist, all.peaks.split, SIMPLIFY=FALSE)
    
    listcode <- strsplit(listcode, "")
    names(listcode) <- listname
    listcode <- sapply(listcode, function(.ele) sum(as.numeric(.ele))==2)
    overlappingPeaks <- peaklist[listcode[names(peaklist)]]
    Peaks <- venn_cnt$Peaks
    names(Peaks) <- gsub(NAME_conn_string, NAME_short_string, names(Peaks))
    overlappingPeaks <- lapply(overlappingPeaks, function(.ele){
        peakNames <- unlist(.ele$peakNames)
        ps <- Peaks[peakNames]
        if(!is.null(ps$old_strand_HH)){
            strand(ps) <- ps$old_strand_HH
            ps$old_strand_HH <- NULL
        }
        gp <- gsub(paste(NAME_short_string, ".*$", sep=""), "", peakNames)
        ps <- split(ps, gp)
        ol <- findOverlaps(query=ps[[1]], subject=ps[[2]], 
                           maxgap=maxgap, minoverlap=minoverlap, 
                           ignore.strand=ignore.strand)
        q <- ps[[1]][queryHits(ol)]
        s <- ps[[2]][subjectHits(ol)]
        cl <- getRelationship(q, s)[,c("insideFeature", "shortestDistance")]
        colnames(cl)[grepl("insideFeature", colnames(cl))] <- "overlapFeature"
        correlation <- cbind(peaks1=names(q), as.data.frame(unname(q)), 
                             peaks2=names(s), as.data.frame(unname(s)), 
                             cl)
        rownames(correlation) <- make.names(paste(names(q), names(s), sep="_"))
        correlation <- correlation[correlation[,"shortestDistance"]<maxgap | 
                                       correlation[,"overlapFeature"] %in% 
                                       c("includeFeature", "inside",
                                         "overlapEnd", "overlapStart"),]
    })
    PeaksList <- sapply(PeaksList, trimPeakList, by="region",
           ignore.strand=ignore.strand,
           keepMetadata=TRUE)
    for(i in 1:n){
        names(PeaksList[[i]]) <- 
            paste(names[i], names(PeaksList[[i]]), sep=NAME_short_string)
    }
    names(PeaksList) <- NULL
    structure(list(venn_cnt=venn_cnt$venn_cnt, 
                   peaklist=peaklist, 
                   overlappingPeaks=overlappingPeaks,
                   all.peaks=PeaksList), 
              class="overlappingPeaks")
}