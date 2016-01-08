annoPeaks <- function(peaks, annoData, 
                      bindingType=c("bothSidesNSS",
                                    "startSite", "endSite", "fullRange",
                                    "bothSidesNearest"), 
                      bindingRegion=c(-5000, 500), 
                      ignore.peak.strand=TRUE,
                      ...){
    bindingType <- match.arg(bindingType)
    stopifnot(inherits(peaks, "GRanges"))
    stopifnot(inherits(annoData, c("annoGR", "GRanges")))
    stopifnot(seqlevelsStyle(peaks)==seqlevelsStyle(annoData))
    stopifnot(length(bindingRegion)==2)
    stopifnot(bindingRegion[1]<=0 && bindingRegion[2]>=1)
    if(ignore.peak.strand){
        peaks$peakstrand <- strand(peaks)
        strand(peaks) <- "*"
    }
    tmp <- annoData
    annotation <- switch(bindingType,
                         startSite={
                             idx <- as.character(strand(tmp))=="-"
                             start(tmp)[idx] <- end(tmp)[idx]
                             width(tmp) <- 1
                             tmp
                         },
                         endSite={
                             idx <- as.character(strand(tmp))!="-"
                             start(tmp)[idx] <- end(tmp)[idx]
                             width(tmp) <- 1
                             tmp
                         },
                         fullRange=annoData,
                         bothSidesNearest=annoData,
                         bothSidesNSS=annoData,
                         annoData)
    rm(tmp)
    if(bindingType %in% c("bothSidesNearest", "bothSidesNSS")){
        extGR <- function(a, b){
            str_a <- as.character(strand(a))!="-"
            s1 <- ifelse(str_a, start(a)+b[1], start(a)-b[2])
            s1[s1<1] <- 1
            start(a) <- s1
            end(a) <- ifelse(str_a, end(a)+b[2], end(a)-b[1])
            a
        }
        peaks.tmp <- extGR(peaks, bindingRegion)
        ol <- findOverlaps(query=peaks.tmp, subject=annotation,
                           maxgap=0L, minoverlap=1L,
                           type="any", select="all",
                           ignore.strand=FALSE)
    }else{
        idx <- as.character(strand(annotation))!="-"
        s1 <- ifelse(idx, start(annotation)+bindingRegion[1], 
                                    start(annotation)+1-bindingRegion[2])
        s1[s1<1] <- 1
        start(annotation) <- s1
        e1 <- ifelse(idx, end(annotation)-1+bindingRegion[2],
                                  end(annotation)-bindingRegion[1])
        seql <- seqlengths(annotation)
        if(length(seql)>0){
            e1_seql <- seql[as.character(seqnames(annotation))]
            id <- e1_seql<e1
            id <- id[!is.na(id)]
            if(length(id)>0) e1[id] <- e1_seql[id]
        }
        end(annotation) <- e1
        if(bindingType=="startSite"){##make sure the downstream is inside gene
            start(annotation)[!idx] <- 
                ifelse(start(annotation)[!idx]>start(annoData)[!idx], 
                       start(annotation)[!idx], 
                       start(annoData)[!idx])
            end(annotation)[idx] <- 
                ifelse(end(annotation)[idx]<end(annoData)[idx], 
                       end(annotation)[idx], 
                       end(annoData)[idx])
        }else{
            if(bindingType=="endSite"){
                start(annotation)[idx] <- 
                    ifelse(start(annotation)[idx]>start(annoData)[idx], 
                           start(annotation)[idx], 
                           start(annoData)[idx])
                end(annotation)[!idx] <- 
                    ifelse(end(annotation)[!idx]<end(annoData)[!idx], 
                           end(annotation)[!idx], 
                           end(annoData)[!idx])
            }
        }
        
        ol <- findOverlaps(query=peaks, subject=annotation,
                           maxgap=0L, minoverlap=1L,
                           type="any", select="all",
                           ignore.strand=FALSE)
    }
    if(length(ol)<1){
        return(GRanges())
    }
    peaks <- peaks[queryHits(ol)]
    anno <- annoData[subjectHits(ol)]
    if(bindingType %in% c("bothSidesNearest", "bothSidesNSS")){
        relations <- getRelationship(peaks, anno)
        ##filter resutls and save the nearest
        keep <- rep(FALSE, length(peaks))
        anno.strand <- as.character(strand(anno))!="-"
        if(bindingType=="bothSidesNSS"){
            keep[relations$insideFeature %in% 
                     c("includeFeature", "overlap")] <- TRUE
            keep.left <- (relations$insideFeature %in% 
                c("upstream", "overlapStart") & !anno.strand) |
                (relations$insideFeature %in% "inside") | 
                ((relations$insideFeature %in% "overlapEnd") & anno.strand)
            keep.right <- (relations$insideFeature %in% 
                c("upstream", "overlapStart") & anno.strand) |
                (relations$insideFeature %in% "inside") |
                ((relations$insideFeature %in% "overlapEnd") & !anno.strand)
            shortestDist <- relations$distanceToStart
        }else{
            keep[relations$insideFeature %in% 
                     c("includeFeature", "inside", "overlap",
                       "overlapEnd", "overlapStart")] <- TRUE
            keep.left <- (relations$insideFeature=="downstream" & anno.strand) |
                (relations$insideFeature=="upstream" & !anno.strand)
            keep.right <- (relations$insideFeature=="upstream" & anno.strand) |
                (relations$insideFeature=="downstream" & !anno.strand)
            shortestDist <- relations$shortestDistance
        }
        names(shortestDist) <- 1:length(peaks)
        whichismin <- function(.ele){
            as.numeric(names(.ele)[.ele==min(.ele)])
        }
        if(sum(keep.left)>=1){
            nearest.left <- tapply(shortestDist[keep.left], 
                                       queryHits(ol)[keep.left],
                                       whichismin, simplify=FALSE)
            keep[unlist(nearest.left)] <- TRUE
        }
        if(sum(keep.right)>=1){
            nearest.right <- 
                tapply(shortestDist[keep.right], 
                       queryHits(ol)[keep.right],
                       whichismin, simplify=FALSE)
            keep[unlist(nearest.right)] <- TRUE
        }
        peaks <- peaks[keep]
        anno <- anno[keep]
    }
    if(!is.null(names(anno))){
        peaks$feature <- names(anno)
    }
    peaks$feature.ranges <- unname(ranges(anno))
    peaks$feature.strand <- strand(anno)
    peaks$distance <- distance(peaks, anno, ignore.strand=FALSE)
    relations <- getRelationship(peaks, anno)
    #peaks$binding.site <- relations$insideFeature
    peaks$insideFeature <- relations$insideFeature
    peaks$distanceToStart <- relations$distanceToStart
    if(ignore.peak.strand){
        strand(peaks) <- peaks$peakstrand
        peaks$peakstrand <- NULL
    }
    mcols(peaks) <- cbind(mcols(peaks), mcols(anno))
    peaks
}