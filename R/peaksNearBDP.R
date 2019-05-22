peaksNearBDP <- function(myPeakList, AnnotationData,
                         MaxDistance=5000L, ...){
    if (missing(myPeakList)) {
        stop("Missing required argument myPeakList!")
    }
    if (!inherits(myPeakList, c("RangedData", "GRanges"))) {
        stop("myPeakList needs to be GRanges object")
        if(is(myPeakList, "RangedData"))
            myPeakList <- toGRanges(myPeakList, format="RangedData")
    }
    if (!missing(AnnotationData)){        
        if (!inherits(AnnotationData, c("RangedData", "GRanges", "annoGR"))) {
            stop("AnnotationData needs to be GRanges or annoGR object")
        }
        if(class(AnnotationData)=="annoGR")
            AnnotationData <- AnnotationData@gr
        if(is(AnnotationData, "RangedData")) 
            AnnotationData <- toGRanges(AnnotationData, format="RangedData")
    }else{
        stop("Missing required argument AnnotationData!")
    }
    stopifnot(length(intersect(seqlevelsStyle(myPeakList),seqlevelsStyle(AnnotationData)))>0)
    stopifnot(is.numeric(MaxDistance))
    
    MaxDistance <- round(MaxDistance[1])
    myPeakList <- unique(myPeakList)
    myPeakList$bdp_idx <- 1:length(myPeakList)
    anno <- annoPeaks(myPeakList, AnnotationData, 
                      bindingType = "nearestBiDirectionalPromoters",
                      bindingRegion = c(-1*MaxDistance, MaxDistance))
    if(length(anno)<1){
        return(list(peaksWithBDP=anno,
                    percentPeaksWithBDP=0,
                    n.peaks=length(myPeakList),
                    n.peaksWithBDP=0))
    }
    anno.s <- split(anno, anno$bdp_idx)
    len <- elementNROWS(anno.s)
    anno.s <- anno.s[len>=2]
    len <- sapply(anno.s, function(.ele){
        std <- .ele$feature.strand
        all(c("+", "-") %in% as.character(.ele$feature.strand))
    })
    anno.s <- anno.s[len]
    list(peaksWithBDP=anno.s,
         percentPeaksWithBDP = length(anno.s)/length(myPeakList),
         n.peaks=length(myPeakList),
         n.peaksWithBDP=length(anno.s))
}
