formatStrand <- function(strand){
  strand <- as.character(strand)
  strand.levels=levels(as.factor(strand))
  strand.allowed.characters=c("1","-1","+","-", "*")
  if(any(!strand.levels %in% strand.allowed.characters))
  {
    warning("All the characters for strand, 
            other than '1', '-1', '+', '-' and '*', 
            will be converted into '*'.")
  }
  strand[strand== "-1"] <- "-"
  strand[strand== "1"] <- "+"
  strand[!strand %in% strand.allowed.characters] <- "*"
  strand
}
###clear seqnames, the format should be chr+NUM
formatSeqnames <- function(gr) {
    if(class(gr)=="GRanges"){
        seqlevels(gr)[grepl("^(\\d+|V?I{0,3}|IV|MT|M|X|Y)$", seqlevels(gr))] <-
            paste("chr", 
                  seqlevels(gr)[grepl("^(\\d+|V?I{0,3}|IV|MT|M|X|Y)$", 
                                      seqlevels(gr))], sep="")
        seqlevels(gr)[seqlevels(gr)=="chrMT"] <- "chrM" 
    }else{
        seqnames(gr)[grepl("^(\\d+|V?I{0,3}|IV|MT|M|X|Y)$", seqnames(gr))] <-
            paste("chr", 
                  seqnames(gr)[grepl("^(\\d+|V?I{0,3}|IV|MT|M|X|Y)$", 
                                     seqnames(gr))], sep="")
        seqnames(gr)[seqnames(gr)=="chrMT"] <- "chrM" 
    }
#    if(seqlevelsStyle(gr)!="UCSC") seqlevelsStyle(gr) <- "UCSC"
    gr
}
RangedData2GRanges <- function(rd){
  ss <- grep("^strand$", colnames(rd), ignore.case=TRUE)
  if(length(ss)>1) 
      stop("input RangedData has multiple columns for strand information")
  if(length(ss)==1){
    rd.strand <- formatStrand(rd[[ss]])
  }else{
    rd.strand <- "*"
  }
  for(coln in colnames(rd)){##because the bug of as(rd, "GRanges")
    rd[[coln]] <- NULL
  }
  rd$strand <- rd.strand
  gr <- as(rd, "GRanges")
  gr
}
getRelationship <- function(queryHits, subjectHits){
    if(!inherits(queryHits, "GRanges")) 
        stop("queryHits must be an object of GRanges")
    if(!inherits(subjectHits, "GRanges")) 
        stop("subjectHits must be an object of GRanges")
    strand <- strand(subjectHits)=="-"
    FeatureStart <- as.numeric(ifelse(strand, 
                                      end(subjectHits), 
                                      start(subjectHits)))
    FeatureEnd <- as.numeric(ifelse(strand, 
                                    start(subjectHits), 
                                    end(subjectHits)))
    PeakStart <- as.numeric(ifelse(strand, end(queryHits), 
                                   start(queryHits)))
    PeakEnd <- as.numeric(ifelse(strand, start(queryHits), end(queryHits)))
    ss <- PeakStart - FeatureStart
    ee <- PeakEnd - FeatureEnd
    se <- PeakStart - FeatureEnd
    es <- PeakEnd - FeatureStart
    upstream <- ifelse(strand, es>0, es<0)
    downstream <- ifelse(strand, se<0, se>0)
    includeFeature <- ifelse(strand, (ss>=0 & ee<=0), (ss<=0 & ee>=0))
    overlap <- ss==0 & ee==0
    inside <- ifelse(strand, (ss<=0 & ee>=0), (ss>=0 & ee<=0))
    overlapStart <- 
        ifelse(strand, (ss>0 & es<=0 & ee>=0), (ss<0 & es>=0 & ee<=0))
    overlapEnd <- 
        ifelse(strand, (ss<0 & se>=0 & ee<=0), (ss>0 & se<=0 & ee>=0))
    insideFeature <- rep(NA, length(queryHits))
    insideFeature[as.logical(includeFeature)] <- "includeFeature"
    insideFeature[as.logical(inside)] <- "inside"
    insideFeature[as.logical(overlap)] <- "overlap"
    insideFeature[as.logical(overlapEnd)] <- "overlapEnd"
    insideFeature[as.logical(downstream)] <- "downstream"
    insideFeature[as.logical(overlapStart)] <- "overlapStart"
    insideFeature[as.logical(upstream)] <- "upstream"
    shortestDistance <- apply(cbind(ss, ee, se, es), 1,
                              function(.ele) min(abs(.ele)))
    shortestDistanceToStart <- apply(cbind(ss, es), 1, 
                                     function(.ele) min(abs(.ele)))
    data.frame(insideFeature=insideFeature, 
               shortestDistance=shortestDistance, 
               ss=ss,
               distanceToStart=shortestDistanceToStart)
}

vennCounts <- function(PeaksList, n, names,
                       maxgap=0L, minoverlap=1L, 
                       by=c("region", "feature", "base"), 
                       ignore.strand=TRUE, 
                       connectedPeaks=c("min", "merge", "keepAll")){
    NAME_conn_string <- "___conn___"
    PeaksList<-lapply(PeaksList, function(Peaks){
        if (inherits(Peaks, "RangedData"))
            Peaks <- toGRanges(Peaks, format="RangedData")
        if (!inherits(Peaks, "GRanges")) {
            stop("No valid Peaks passed in. It needs to be GRanges object")
        }
        if(ignore.strand) {
            .gr <- paste(seqnames(Peaks), start(Peaks), end(Peaks))
        }else{
            .gr <- paste(seqnames(Peaks), start(Peaks), 
                         end(Peaks), strand(Peaks))
        }
        if(any(duplicated(.gr)))
            stop("Inputs contains duplicated ranges. 
                 please recheck your inputs.")
        if(any(is.null(names(Peaks))) || 
               any(is.na(names(Peaks))) || 
               any(duplicated(names(Peaks)))) {
            message("duplicated or NA names found. 
                    Rename all the names by numbers.")
            names(Peaks) <- formatC(1:length(Peaks), 
                                    width=nchar(length(Peaks)), 
                                    flag='0')
        }
        feature <- mcols(Peaks)$feature
        mcols(Peaks) <- NULL
        if(by=="feature") {
            if(is.null(feature)) 
                stop("Need feature metadata for each inputs")
            mcols(Peaks)$feature <- feature
        }
        Peaks
    })
    if(by=="base"){
        # try coverage and then split
        # problem for coverage: can not seperate for different source
        for(i in 1:n){
            names(PeaksList[[i]]) <- 
                paste(names[i], names(PeaksList[[i]]), sep=NAME_conn_string)
        }
        names(PeaksList) <- NULL
        Peaks <- unlist(GRangesList(PeaksList))
        #if(ignore.strand) strand(Peaks) <- "*"
        #cov <- coverage(Peaks)
        disj <- disjoin(Peaks, ignore.strand=ignore.strand)
        ol <- findOverlaps(Peaks, disj, ignore.strand=ignore.strand)
        ol.query <- Peaks[queryHits(ol)]
        ol.subject <- disj[subjectHits(ol)]
        group <- gsub(paste(NAME_conn_string, ".*$", sep=""), 
                      "", 
                      names(ol.query))
        subject <- subjectHits(ol)
        ncontrasts <- n
        noutcomes <- 2^ncontrasts
        outcomes <- matrix(0,noutcomes,ncontrasts)
        colnames(outcomes) <- names
        for (j in 1:ncontrasts)
            outcomes[,j] <- rep(0:1,times=2^(j-1),each=2^(ncontrasts-j))
        gps <- split(group, subject)
        xlist <- list()
        for(i in 1:ncontrasts)
            xlist[[i]] <- 
            factor(as.numeric(sapply(gps, 
                                     function(.ele) 
                                         names[ncontrasts-i+1] %in% .ele)), 
                   levels=c(0,1))
        counts <- do.call(cbind, xlist)
        counts <- counts[, ncontrasts:1]
        counts <- counts - 1
        counts <- apply(counts, 1, paste, collapse="")
        idx <- apply(outcomes, 1, paste, collapse="")
        wids <- width(disj[as.numeric(names(gps))])
        wids <- split(wids, counts)
        wids <- wids[idx]
        names(wids) <- idx
        counts <- sapply(wids, sum)
        venn_cnt <- structure(cbind(outcomes, Counts=counts), 
                              class="VennCounts")
        return(list(venn_cnt=venn_cnt, xlist=NULL, 
                    PeaksList=PeaksList, all=all, Peaks=Peaks))
    }
    if(by=="feature"){
        features <- lapply(PeaksList, 
                           function(.ele) unique(as.character(.ele$feature)))
        all_features <- unique(unlist(features))
        all <- lapply(features, function(.f)  all_features %in% .f)
        all <- do.call(cbind, all)
        rownames(all) <- all_features
        for(j in 1:n){
            all[,j] <- ifelse(all[,j], names[j], NA)
        }
        all <- split(all, rownames(all))
        Peaks <- NULL
    }else{
        ##get all merged peaks
        for(i in 1:n){
            names(PeaksList[[i]]) <- 
                paste(names[i], names(PeaksList[[i]]), sep=NAME_conn_string)
        }
        names(PeaksList) <- NULL
        Peaks <- unlist(GRangesList(PeaksList))
        if(ignore.strand) {
            Peaks$old_strand_HH <- strand(Peaks)
            strand(Peaks) <- "*"
        }
        
        if(length(Peaks)<10000){
            ol <- as.data.frame(findOverlaps(Peaks, maxgap=maxgap, 
                                             minoverlap=minoverlap, 
                                             select="all",
                                             ignoreSelf=TRUE, 
                                             ignoreRedundant=TRUE))
            ##all connected peaks
            olm <- cbind(names(Peaks[ol[,1]]), names(Peaks[ol[,2]]))
            edgeL <- c(split(olm[,2], olm[,1]), split(olm[,1], olm[,2]))
            nodes <- unique(as.character(olm))
            ##use graph to extract all the connected peaks
            gR <- new("graphNEL", nodes=nodes, edgeL=edgeL)
            Merged <- connectedComp(ugraph(gR))
        }else{
            Peaks.list <- split(Peaks, seqnames(Peaks))
            Merged <- lapply(Peaks.list, function(.peaks.list){
                .ol <- findOverlaps(.peaks.list,
                                    maxgap=maxgap, minoverlap=minoverlap, 
                                    select="all",
                                    ignoreSelf=TRUE, ignoreRedundant=TRUE)
                olm <- cbind(names(.peaks.list[queryHits(.ol)]), 
                             names(.peaks.list[subjectHits(.ol)]))
                edgeL <- c(split(olm[,2], olm[,1]), split(olm[,1], olm[,2]))
                nodes <- unique(as.character(olm))
                ##use graph to extract all the connected peaks
                gR <- new("graphNEL", nodes=nodes, edgeL=edgeL)
                connectedComp(ugraph(gR))
            })
            Merged <- unlist(Merged, recursive=FALSE)
            nodes <- unique(as.character(unlist(Merged)))
            rm(Peaks.list)
        }
        
        Left <- as.list(names(Peaks)[!names(Peaks) %in% nodes])
        all <- c(Merged, Left)
    }
    ##venn count
    ncontrasts <- n
    noutcomes <- 2^ncontrasts
    outcomes <- matrix(0,noutcomes,ncontrasts)
    colnames(outcomes) <- names
    for (j in 1:ncontrasts)
        outcomes[,j] <- rep(0:1,times=2^(j-1),each=2^(ncontrasts-j))
    
    xlist <- list()
    xlist1 <- list()
    ## time consuming step, FIXME!!
    ## Fixed by paste connection string before lapply 
    NAME_conn_string_wild <- paste(NAME_conn_string, ".*?$", sep="")
    for (i in 1:ncontrasts){
        NAME_conn_string_contrasts_wild <- paste("^", 
                                                 names[ncontrasts-i+1], 
                                                 NAME_conn_string,
                                                 sep="")
        xlist[[i]] <- factor(as.numeric(unlist(lapply(all, function(.ele) 
            any(grepl(NAME_conn_string_contrasts_wild, .ele))))),
            levels=c(0,1))
        if(connectedPeaks=="merge"){
            xlist1[[i]] <- xlist[[i]]
        }else{
            xlist1[[i]] <- 
                factor(as.numeric(unlist(lapply(all, function(.ele) {
                ##count involved nodes in each group
                .ele <- gsub(NAME_conn_string_wild, "", .ele)
                .ele <- table(.ele)
                rep(names[ncontrasts-i+1] %in% names(.ele), min(.ele))
            }))), levels=c(0,1))
        }
    }
    
    counts <- as.vector(table(xlist1))
    venn_cnt <- structure(cbind(outcomes, Counts=counts), class="VennCounts")
    
    if(connectedPeaks=="keepAll"){
        NAME_conn_string_wild <- paste(NAME_conn_string, ".*$", sep="")
        all.m <- lapply(all, 
                        function(.ele){gsub(NAME_conn_string_wild, "", .ele)})
        all.m <- all.m[sapply(all.m, function(.ele) length(unique(.ele))>1)]
        all.count <- ifelse(rowSums(outcomes)>1, NA, counts)
        all.count <- all.count * outcomes
        for(j in 1:ncol(all.count)){
            ylist <- list()
            for(i in 1:ncontrasts){
                ylist[[i]] <- 
                    factor(as.numeric(unlist(lapply(all.m, function(.ele){
                    .ele <- table(.ele)
                    times <- .ele[colnames(all.count)[j]]
                    if(is.na(times)) times <- 0
                    rep(names[ncontrasts-i+1] %in% names(.ele), times)
                }))), levels=c(0,1))
            }
            all.count[is.na(all.count[,j]),j] <- 
                as.vector(table(ylist))[is.na(all.count[,j])]
        }
        colnames(all.count) <- paste("count", colnames(all.count), sep=".")
        venn_cnt <- structure(cbind(outcomes, Counts=counts, all.count), 
                              class="VennCounts")
    }
    
    return(list(venn_cnt=venn_cnt, xlist=xlist, 
                PeaksList=PeaksList, all=all, 
                Peaks=Peaks))
}