##require library graph, RBGL
findOverlappingPeaks <- function(Peaks1, Peaks2, maxgap = -1L,minoverlap=0L,
                                 multiple = c(TRUE, FALSE),
                                 NameOfPeaks1 = "TF1", NameOfPeaks2 = "TF2", 
                                 select=c("all", "first", "last", "arbitrary"), 
                                 annotate =0, ignore.strand=TRUE, 
                                 connectedPeaks=c("min", "merge"), ...){
        .Deprecated("findOverlapsOfPeaks")
        ###check inputs
        NAME_conn_string <- "___conn___"
        NAME_short_string <- "__"
        NAME_long_string <- "///"
        PeaksList <- list(...)
        PeaksList<-lapply(PeaksList, function(Peaks){
            if (inherits(Peaks, "RangedData"))
                Peaks <- toGRanges(Peaks, format="RangedData")
            if (!inherits(Peaks, "GRanges")) {
                stop("No valid Peaks passed in. It needs to be GRanges object")
            }
            if(any(is.na(names(Peaks))) || any(duplicated(names(Peaks)))) {
                message("duplicated or NA names found. 
                        Rename all the names by numbers.")
                names(Peaks) <- formatC(1:nrow(data), 
                                        width=nchar(nrow(data)), 
                                        flag='0')
            }
            Peaks
        })
        n <- length(PeaksList)
        if(n>0){
            if(n==1){
                PeaksList <- PeaksList[[1]]
                n <- length(PeaksList)
                names <- names(PeaksList)
                if(is.null(names)) names <- paste("peaks", 1:n, sep="")
            }else{
                ##save dots arguments names
                dots <- substitute(list(...))[-1]
                names <- unlist(sapply(dots, deparse))
            }
        }else{
            names <- NULL
        }
        flagRD <- FALSE
        if((!missing(Peaks1)) || (!missing(Peaks2))){
            if(!missing(Peaks2)){##
                if(inherits(Peaks2, "RangedData")){
                    Peaks2 <- toGRanges(Peaks2, format="RangedData")
                    flagRD <- TRUE
                }
                if (!inherits(Peaks2, "GRanges")) {
                    stop("No valid Peaks passed in. 
                         It needs to be GRanges object")
                }
                if(n!=0) PeaksList <- c(list(Peaks2), PeaksList)
                else PeaksList <- list(Peaks2)
                n <- n+1
                if(missing(NameOfPeaks2)) NameOfPeaks2 <- "Peaks2"
                names <- c(NameOfPeaks2, names)
            }
            if(!missing(Peaks1)){
                if(inherits(Peaks1, "RangedData")){
                    Peaks1 <- toGRanges(Peaks1, format="RangedData")
                    flagRD <- TRUE
                }
                if (!inherits(Peaks1, "GRanges")) {
                    stop("No valid Peaks passed in. 
                         It needs to be GRanges object")
                }
                if(n!=0) PeaksList <- c(list(Peaks1), PeaksList)
                else PeaksList <- list(Peaks1)
                n <- n+1
                if(missing(NameOfPeaks1)) NameOfPeaks1 <- "Peaks1"
                names <- c(NameOfPeaks1, names)
            }
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
        ##handle colnames of metadata
        metacolnames <- lapply(PeaksList, function(Peaks)
            colnames(mcols(Peaks)))
        metacolnames <- Reduce(intersect, metacolnames)
        metacolclass <- do.call(rbind, lapply(PeaksList, function(Peaks)
            sapply(mcols(Peaks)[, metacolnames, drop=FALSE], class)))
        metacolclass <- apply(metacolclass, 2, 
                              function(.ele) length(unique(.ele))==1)
        metacolnames <- metacolnames[metacolclass]
        PeaksList <- lapply(PeaksList, function(Peaks){
            mcols(Peaks) <- mcols(Peaks)[, metacolnames]
            Peaks
        })
        ##get all merged peaks
        for(i in 1:n){
            names(PeaksList[[i]]) <- 
                paste(names[i], names(PeaksList[[i]]), 
                      sep=NAME_conn_string)
        }
        Peaks <- unlist(GRangesList(PeaksList))
        if(ignore.strand) strand(Peaks) <- "*"
        
        ol <- as.data.frame(findOverlaps(Peaks, maxgap=maxgap, 
                                         minoverlap=minoverlap, 
                                         select="all",
                                         drop.self=TRUE, 
                                         drop.redundant=TRUE))
        olm <- cbind(names(Peaks[ol[,1]]), names(Peaks[ol[,2]]))
        edgeL <- c(split(olm[,2], olm[,1]), split(olm[,1], olm[,2]))
        nodes <- unique(as.character(olm))
        ##use graph to extract all the connected peaks
        gR <- new("graphNEL", nodes=nodes, edgeL=edgeL)
        Merged <- connectedComp(ugraph(gR))        
        Left <- as.list(names(Peaks)[!names(Peaks) %in% nodes])
        all <- c(Merged, Left)
        
        ##venn count
        ncontrasts <- n
        noutcomes <- 2^ncontrasts
        outcomes <- matrix(0,noutcomes,ncontrasts)
        colnames(outcomes) <- names
        for (j in 1:ncontrasts)
            outcomes[,j] <- rep(0:1,times=2^(j-1),each=2^(ncontrasts-j))
        xlist <- list()
        xlist1 <- list()
        for (i in 1:ncontrasts){
            xlist[[i]] <- factor(as.numeric(unlist(lapply(all, function(.ele) 
                any(grepl(paste("^",
                                names[ncontrasts-i+1], 
                                NAME_conn_string, sep=""), .ele))))),
                                 levels=c(0,1))
            if(connectedPeaks=="merge"){
                xlist1[[i]] <- xlist[[i]]
            }else{
                xlist1[[i]] <- 
                    factor(as.numeric(unlist(lapply(all, function(.ele) {
                    ##count involved nodes in each group
                    if(length(.ele)>2){
                        .ele <- gsub(paste(NAME_conn_string, ".*?$", sep=""), 
                                     "", .ele)
                        .ele <- table(.ele)
                        rep(names[ncontrasts-i+1] %in% names(.ele), 
                            min(.ele))
                    }else{
                        any(grepl(paste("^",names[ncontrasts-i+1], 
                                        NAME_conn_string, sep=""), 
                                  .ele))
                    }
                }))), levels=c(0,1))
            }
        }
        counts <- as.vector(table(xlist1))
        venn_cnt <- structure(cbind(outcomes, Counts=counts), 
                              class="VennCounts")
        xlist <- do.call(rbind, xlist)
        xlist <- xlist - 1
        xlist <- xlist[nrow(xlist):1,,drop=FALSE] 
        ## reverse xlist to match the order of names
        xlist <- apply(xlist, 2, paste, collapse="")
        all <- do.call(rbind, mapply(function(.ele, .id) cbind(.id, .ele), 
                                     all, 1:length(all), SIMPLIFY=FALSE))
        all.peaks <- Peaks[all[,2]]
        all.peaks$gpForFindOverlapsOfPeaks <- all[, 1]
        all.peaks.rd <- reduce(all.peaks, min.gapwidth=maxgap+1L, with.revmap=TRUE)
        mapping <- all.peaks.rd$revmap
        m <- sapply(mapping, length)
        mIndex <- rep(1:length(mapping), m)
        mLists <- unlist(mapping)
        mcols <- mcols(all.peaks[mLists])
        mcols$peakNames <- gsub(NAME_conn_string, 
                                NAME_short_string, 
                                names(all.peaks[mLists]))
        mcolsn <- sapply(mcols[1, ], class)
        mapping <- DataFrame(HHH_row___H=1:length(all.peaks.rd))
        for(.name in names(mcolsn)){
            .dat <- split(mcols[, .name], mIndex)
            mapping[, .name] <- switch(mcolsn[.name],
                                       logical=LogicalList(.dat),
                                       integer=IntegerList(.dat),
                                       numeric=NumericList(.dat),
                                       character=CharacterList(.dat),
                                       rle=RleList(.dat),
                                       ComplexList(.dat))
        }
        mapping$HHH_row___H <- NULL
        mcols(all.peaks.rd) <- mapping
        names(all.peaks.rd) <- sapply(all.peaks.rd$peakNames, 
                                      paste, 
                                      collapse=NAME_short_string)
        all.peaks.rd$gpForFindOverlapsOfPeaks <- 
            unlist(lapply(all.peaks.rd$gpForFindOverlapsOfPeaks, unique)) 
        all <- split(all.peaks.rd, all.peaks.rd$gpForFindOverlapsOfPeaks)
        all <- all[order(as.numeric(names(all)))] 
        ##important, and length(all)==length(xlist)
        if(length(all)!=length(xlist)) 
            stop("length of all should be equal to length of xlist. 
                 Please report the bug. Thanks.")
        listname <- apply(outcomes, 1, 
                          function(id) paste(names[as.logical(id)], 
                                             collapse=NAME_long_string))
        listcode <- apply(outcomes, 1, paste, collapse="")
        listname <- listname[-1]
        listcode <- listcode[-1]
        peaklist <- list()
        for(i in 1:length(listcode)){
            sublist <- all[xlist==listcode[i]]
            if(length(sublist)>0) 
                peaklist[[listname[i]]]<-unlist(sublist, use.names=FALSE)
        }
        correlation <- list()
        names(Peaks) <- gsub(NAME_conn_string, NAME_short_string, names(Peaks))
        npl <- names(peaklist)
        for(i in 1:length(npl))
        {
            npln <- unlist(strsplit(npl[i], NAME_long_string))
            if(length(npln)==2){
                pl <- peaklist[[i]]$peakNames
                pl <- pl[sapply(pl, length)>1]
                pl1 <- lapply(pl, 
                              function(.ele) 
                                  .ele[grepl(paste("^", npln[1], 
                                                   NAME_short_string, 
                                                   sep=""), 
                                             .ele)][1])
                pl2 <- lapply(pl, 
                              function(.ele) 
                                  .ele[grepl(paste("^", npln[2], 
                                                   NAME_short_string, 
                                                   sep=""), .ele)][1])
                idsel <- mapply(function(.p1, .p2) is.na(.p1)+is.na(.p2)==0, 
                                pl1, pl2)
                pl1 <- Peaks[unlist(pl1[idsel])]
                pl2 <- Peaks[unlist(pl2[idsel])]
                cl <- 
                    getRelationship(pl1, pl2)[,c("insideFeature", 
                                                 "shortestDistance")]
                colnames(cl)[grepl("insideFeature", colnames(cl))] <- 
                    "overlapFeature"
                correlation[[npl[i]]] <- 
                    cbind(peaks1=names(pl1), as.data.frame(pl1), 
                          peaks2=names(pl2), as.data.frame(pl2), 
                          cl)
                rownames(correlation[[npl[i]]]) <- 
                    paste(names(pl1), names(pl2), sep="_")
            }
        }
        for(i in 1:length(peaklist)){
            peaklist[[i]]$gpForFindOverlapsOfPeaks <- NULL
        }
        
        if((!missing(Peaks1)) && (!missing(Peaks2))){
            sampleName <- 
                npl==paste(NameOfPeaks1, 
                           NAME_long_string, 
                           NameOfPeaks2, sep="")
            if(!any(sampleName)) 
                sampleName <- npl==paste(NameOfPeaks2, 
                                         NAME_long_string, 
                                         NameOfPeaks1, sep="")
            sampleName <- npl[sampleName]
            Peaks1withOverlaps <- 
                Peaks[names(Peaks) %in% correlation[[sampleName]]$peaks1]
            names(Peaks1withOverlaps) <- 
                gsub(paste(NameOfPeaks1, NAME_short_string, sep=""), 
                     "", 
                     names(Peaks1withOverlaps))
            Peaks2withOverlaps <- 
                Peaks[names(Peaks) %in% correlation[[sampleName]]$peaks2]
            names(Peaks2withOverlaps) <- 
                gsub(paste(NameOfPeaks2, NAME_short_string, sep=""), 
                     "", 
                     names(Peaks2withOverlaps))
            if(flagRD){
                peaklist <- lapply(peaklist, as, Class="RangedData")
                correlation <- lapply(correlation, as, Class="RangedData")
                Peaks1withOverlaps <- as(Peaks1withOverlaps, "RangedData")
                Peaks2withOverlaps <- as(Peaks2withOverlaps, "RangedData")
            }
            mergedPeaks <- peaklist[[sampleName]]##To fit old version
            mergedPeaks$peakNames <- NULL
            structure(list(venn_cnt=venn_cnt, 
                           peaklist=peaklist, 
                           overlappingPeaks=correlation,
                           OverlappingPeaks=correlation[[sampleName]],
                           MergedPeaks=mergedPeaks,##
                           Peaks1withOverlaps=Peaks1withOverlaps,
                           Peaks2withOverlaps=Peaks2withOverlaps), 
                      class="overlappingPeaks")
        }else{
            structure(list(venn_cnt=venn_cnt, 
                           peaklist=peaklist, 
                           overlappingPeaks=correlation), 
                      class="overlappingPeaks")
        }
    }
