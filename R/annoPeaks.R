#' Annotate peaks
#' 
#' Annotate peaks by annoGR object in the given range.
#' 
#' 
#' @param peaks peak list, \link[GenomicRanges:GRanges-class]{GRanges} object
#' @param annoData annotation data, \link[GenomicRanges:GRanges-class]{GRanges}
#' object
#' @param bindingType Specifying the criteria to associate peaks with
#' annotation. Here is how to use it together with the parameter bindingRegion
#' \itemize{ \item To obtain peaks within 5kb upstream and up to 3kb downstream
#' of TSS within the gene body, set bindingType = "startSite" and bindingRegion
#' = c(-5000, 3000) \item To obtain peaks up to 5kb upstream within the gene
#' body and 3kb downstream of gene/Exon End, set bindingType = "endSite" and
#' bindingRegion = c(-5000, 3000) \item To obtain peaks from 5kb upstream to
#' 3kb downstream of genes/Exons , set bindingType = "fullRange" and
#' bindingRegion = c(-5000, 3000) \item To obtain peaks with nearest
#' bi-directional promoters within 5kb upstream and 3kb downstream of TSS, set
#' bindingType = "nearestBiDirectionalPromoters" and bindingRegion = c(-5000,
#' 3000) } \describe{ \item{startSite}{start position of the feature (strand is
#' considered)} \item{endSite}{end position of the feature (strand is
#' considered)} \item{fullRange}{whole range of the feature}
#' \item{nearestBiDirectionalPromoters}{nearest promoters from both direction
#' of the peaks (strand is considered). It will report bidirectional promoters
#' if there are promoters in both directions in the given region (defined by
#' bindingRegion). Otherwise, it will report the closest promoter in one
#' direction.} }
#' @param bindingRegion Annotation range used together with bindingType, which
#' is a vector with two integer values, default to c (-5000, 5000). The first
#' one must be no bigger than 0, which means upstream. And the sec ond one must
#' be no less than 1, which means downstream (1 is the site position, 2 is the
#' next base of the site position). For details, see bindingType.
#' @param ignore.peak.strand ignore the peaks strand or not.
#' @param select "all" or "bestOne". Return the annotation containing all or
#' the best one.  The "bestOne" is selected by the shortest distance to the
#' sites and then similarity between peak and annotations.  Ignored if
#' bindingType is nearestBiDirectionalPromoters.
#' @param ...  Not used.
#' @return Output is a GRanges object of the annotated peaks.
#' @export
#' @importFrom GenomeInfoDb seqlevelsStyle seqlengths
#' @importFrom BiocGenerics strand start end width pos
#' @importFrom S4Vectors queryHits subjectHits
#' @author Jianhong Ou
#' @seealso See Also as \code{\link{annotatePeakInBatch}}
#' @keywords misc
#' @examples
#'     library(ensembldb)
#'     library(EnsDb.Hsapiens.v75)
#'     data("myPeakList")
#'     annoGR <- toGRanges(EnsDb.Hsapiens.v75)
#'     seqlevelsStyle(myPeakList) <- seqlevelsStyle(annoGR)
#'     annoPeaks(myPeakList, annoGR)
#' 
annoPeaks <- function(peaks, annoData, 
                      bindingType=c("nearestBiDirectionalPromoters",
                                    "startSite", "endSite", "fullRange"), 
                      bindingRegion=c(-5000, 5000), 
                      ignore.peak.strand=TRUE,
                      select=c("all", "bestOne"), # bestOne will output the one with best score
                      ...){
    select <- match.arg(select)
    if(bindingType[1] %in% 
       c("bothSidesNearest", "nearestBiDirectionalPromoters", "bothSidesNSS")){
        bindingType <- bindingType[1]
        if(bindingType=="bothSidesNSS"){
            bindingType <- "nearestBiDirectionalPromoters"
        }
        if(select!="all"){
            select <- "all"
            message("nearestbiDirectionalPromoters do not support select=bestOne")
        }
    }else{
        bindingType <- match.arg(bindingType)
    }
    stopifnot(inherits(peaks, "GRanges"))
    stopifnot(inherits(annoData, c("annoGR", "GRanges")))
    check_seqlevel <- tryCatch({
    	seqlevelsStyle(peaks)
    	seqlevelsStyle(annoData)
    }, error = function(w) {
    	message("Warning: seqlevel style not recognized, you are probably using custom GRange objects: ", w$message)
    	return(NA)
    })
    if (!any(is.na(check_seqlevel))) {
    	stopifnot(length(intersect(seqlevelsStyle(peaks),seqlevelsStyle(annoData)))>0)
    }
    # stopifnot(length(intersect(seqlevelsStyle(peaks),seqlevelsStyle(annoData)))>0)
    stopifnot(length(bindingRegion)==2)
    stopifnot(bindingRegion[1]<=0 && bindingRegion[2]>=1)
    if(ignore.peak.strand){
        peaks$peakstrand <- strand(peaks)
        strand(peaks) <- "*"
    }
    if(is.null(names(peaks))){
        names(peaks) <- paste0("X", 1:length(peaks))
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
                         nearestBiDirectionalPromoters=annoData,
                         annoData)
    annotation.bck <- annotation
    rm(tmp)
    if(bindingType %in% c("bothSidesNearest", "nearestBiDirectionalPromoters")){
        if(bindingType=="bothSidesNearest"){
            extGR <- function(a, b){
                str_a <- as.character(strand(a))!="-"
                s1 <- ifelse(str_a, start(a)+b[1], start(a)-b[2])
                s1[s1<1] <- 1
                start(a) <- s1
                s2 <- ifelse(str_a, end(a)+b[2], end(a)-b[1])
                s2.idx <- which(s2 > seqlengths(a)[as.character(seqnames(a))])
                if(length(s2.idx)>0){
                    s2[s2.idx] <- seqlengths(a)[as.character(seqnames(a[s2.idx]))]
                }
                end(a) <- s2
                a
            }
            peaks.tmp <- extGR(peaks, bindingRegion)
            ol <- findOverlaps(query=peaks.tmp, subject=annotation,
                               type="any", select="all",
                               ignore.strand=FALSE)
        }else{
            ##bindingType=="nearestBiDirectionalPromoters"
            annotation <- promoters(annotation, 
                                    upstream=-1*bindingRegion[1],
                                    downstream=bindingRegion[2])
            ol <- findOverlaps(query=peaks, subject=annotation,
                               type="any", select="all",
                               ignore.strand=FALSE)
        }
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
                           type="any", select="all",
                           ignore.strand=FALSE)
    }
    if(length(ol)<1){
        return(GRanges())
    }
    peaks <- peaks[queryHits(ol)]
    anno <- annoData[subjectHits(ol)]
    annotation.bck.hits <- annotation.bck[subjectHits(ol)]
    if(bindingType %in% c("bothSidesNearest", "nearestBiDirectionalPromoters")){
        relations <- if(bindingType=="bothSidesNearest") getRelationship(peaks, anno) else getRelationship(peaks, promoters(unname(as(anno, "GRanges")), upstream=0, downstream=1))
        ##filter resutls and save the nearest
        keep <- rep(FALSE, length(peaks))
        anno.strand <- as.character(strand(anno))!="-"
        if(bindingType=="nearestBiDirectionalPromoters"){
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
        annotation.bck.hits <- annotation.bck.hits[keep]
    }
    peaks$peak <- names(peaks)
    if(!is.null(names(anno))){
        peaks$feature <- names(anno)
    }
    peaks$feature.ranges <- unname(ranges(anno))
    peaks$feature.strand <- strand(anno)
    peaks$distance <- distance(peaks, anno, ignore.strand=FALSE)
    relations <- getRelationship(peaks, anno)
    #peaks$binding.site <- relations$insideFeature
    peaks$insideFeature <- relations$insideFeature
    peaks$distanceToSite <- distance(peaks, annotation.bck.hits, 
                                     ignore.strand=ignore.peak.strand)
    if(ignore.peak.strand){
        strand(peaks) <- peaks$peakstrand
        peaks$peakstrand <- NULL
    }
    mcols(peaks) <- cbind(mcols(peaks), mcols(anno))
    if(select=="bestOne"){
        if(length(peaks)==0) return(peaks)
        annoscore <- -1 * annoScore(peaks, anno)
        peaks$ANNOPEAKS__peak.oid <- 1:length(peaks)
        peaks <- peaks[order(peaks$peak, peaks$distanceToSite, annoscore)]
        peaks <- peaks[!duplicated(peaks$peak)]
        peaks <- peaks[order(peaks$ANNOPEAKS__peak.oid)]
        peaks$ANNOPEAKS__peak.oid <- NULL
    }
    peaks
}
