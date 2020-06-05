##require library graph, RBGL


#' Find the overlapping peaks for two peak ranges.
#' 
#' Find the overlapping peaks for two input peak ranges.
#' 
#' The new function findOverlapsOfPeaks is recommended.
#' 
#' Efficiently perform overlap queries with an interval tree implemented in
#' IRanges.
#' 
#' @aliases findOverlappingPeaks findOverlappingPeaks-deprecated
#' @param Peaks1 GRanges: See example below.
#' @param Peaks2 GRanges: See example below.
#' @param maxgap,minoverlap Used in the internal call to \code{findOverlaps()}
#' to detect overlaps. See
#' \code{?\link[IRanges:findOverlaps-methods]{findOverlaps}} in the
#' \pkg{IRanges} package for a description of these arguments.
#' @param multiple TRUE or FALSE: TRUE may return multiple overlapping peaks in
#' Peaks2 for one peak in Peaks1; FALSE will return at most one overlapping
#' peaks in Peaks2 for one peak in Peaks1. This parameter is kept for backward
#' compatibility, please use select.
#' @param NameOfPeaks1 Name of the Peaks1, used for generating column name.
#' @param NameOfPeaks2 Name of the Peaks2, used for generating column name.
#' @param select all may return multiple overlapping peaks, first will return
#' the first overlapping peak, last will return the last overlapping peak and
#' arbitrary will return one of the overlapping peaks.
#' @param annotate Include overlapFeature and shortestDistance in the
#' OverlappingPeaks or not.  1 means yes and 0 means no. Default to 0.
#' @param ignore.strand When set to TRUE, the strand information is ignored in
#' the overlap calculations.
#' @param connectedPeaks If multiple peaks involved in overlapping in several
#' groups, set it to "merge" will count it as only 1, while set it to "min"
#' will count it as the minimal involved peaks in any concered groups
#' @param \dots Objects of \link[GenomicRanges:GRanges-class]{GRanges}: See
#' also \code{\link{findOverlapsOfPeaks}}.
#' @return \item{OverlappingPeaks}{a data frame consists of input peaks
#' information with added information: overlapFeature (upstream: peak1 resides
#' upstream of the peak2; downstream: peak1 resides downstream of the peak2;
#' inside: peak1 resides inside the peak2 entirely; overlapStart: peak1
#' overlaps with the start of the peak2; overlapEnd: peak1 overlaps with the
#' end of the peak2; includeFeature: peak1 include the peak2 entirely) and
#' shortestDistance (shortest distance between the overlapping peaks)}
#' \item{MergedPeaks}{GRanges contains merged overlapping peaks}
#' @author Lihua Julie Zhu
#' @seealso findOverlapsOfPeaks, annotatePeakInBatch, makeVennDiagram
#' @references 1.Interval tree algorithm from: Cormen, Thomas H.; Leiserson,
#' Charles E.; Rivest, Ronald L.; Stein, Clifford. Introduction to Algorithms,
#' second edition, MIT Press and McGraw-Hill. ISBN 0-262-53196-8
#' 
#' 2.Zhu L.J. et al. (2010) ChIPpeakAnno: a Bioconductor package to annotate
#' ChIP-seq and ChIP-chip data. BMC Bioinformatics 2010, 11:237
#' doi:10.1186/1471-2105-11-237
#' 
#' 3. Zhu L (2013). Integrative analysis of ChIP-chip and ChIP-seq dataset.  In
#' Lee T and Luk ACS (eds.), Tilling Arrays, volume 1067, chapter 4, pp. -19.
#' Humana Press. http://dx.doi.org/10.1007/978-1-62703-607-8_8
#' @keywords misc
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @importFrom S4Vectors mcols DataFrame
#' @importClassesFrom graph graphNEL
#' @importFrom graph ugraph
#' @importFrom RBGL connectedComp
#' @importFrom utils data
#' @examples
#' 
#'     if (interactive())
#'     {    
#'     peaks1 = 
#'         GRanges(seqnames=c(6,6,6,6,5), 
#'                 IRanges(start=c(1543200,1557200,1563000,1569800,167889600),
#'                         end=c(1555199,1560599,1565199,1573799,167893599),
#'                         names=c("p1","p2","p3","p4","p5")),
#'                 strand=as.integer(1))
#'     peaks2 = 
#'         GRanges(seqnames=c(6,6,6,6,5), 
#'                 IRanges(start=c(1549800,1554400,1565000,1569400,167888600),
#'                         end=c(1550599,1560799,1565399,1571199,167888999),
#'                         names=c("f1","f2","f3","f4","f5")),
#'                 strand=as.integer(1))
#'     t1 =findOverlappingPeaks(peaks1, peaks2, maxgap=1000, 
#'           NameOfPeaks1="TF1", NameOfPeaks2="TF2", select="all", annotate=1) 
#'     r = t1$OverlappingPeaks
#'     pie(table(r$overlapFeature))
#'     as.data.frame(t1$MergedPeaks)
#'     }
#' 
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
                if(is.null(names)) names <- paste("peaks", seq.int(n), sep="")
            }else{
                ##save dots arguments names
                dots <- substitute(list(...))[-1]
                names <- unlist(sapply(dots, deparse))
            }
        }else{
            names <- NULL
        }
        if((!missing(Peaks1)) || (!missing(Peaks2))){
            if(!missing(Peaks2)){##
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
            sapply(mcols(Peaks)[, metacolnames, drop=FALSE], 
                   function(.ele) class(.ele)[1])))
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
            outcomes[,j] <- rep(0:1,times=2^(j-1), 
                                each=2^(ncontrasts-j))
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
        xlist <- apply(xlist, 2, base::paste, collapse="")
        all <- do.call(rbind, mapply(function(.ele, .id) cbind(.id, .ele), 
                                     all, 1:length(all), SIMPLIFY=FALSE))
        all.peaks <- Peaks[all[,2]]
        all.peaks$gpForFindOverlapsOfPeaks <- all[, 1]
        all.peaks.rd <- reduce(all.peaks, min.gapwidth=maxgap+1L,
                               with.revmap=TRUE)
        mapping <- all.peaks.rd$revmap
        m <- sapply(mapping, length)
        mIndex <- rep(1:length(mapping), m)
        mLists <- unlist(mapping)
        mcols <- mcols(all.peaks[mLists])
        mcols$peakNames <- gsub(NAME_conn_string, 
                                NAME_short_string, 
                                names(all.peaks[mLists]))
        mcolsn <- sapply(mcols[1, ], function(.ele) class(.ele)[1])
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
                                      base::paste,
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
        listcode <- apply(outcomes, 1, base::paste, collapse="")
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
        for(i in seq_along(npl))
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
