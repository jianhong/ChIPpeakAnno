bdp <- function(peaks, annoData, maxgap=2000L, ...){
    stopifnot(inherits(peaks, "GRanges"))
    stopifnot(inherits(annoData, c("annoGR", "GRanges")))
    stopifnot(seqlevelsStyle(peaks)==seqlevelsStyle(annoData))
    stopifnot(is.numeric(maxgap))
    maxgap <- round(maxgap[1])
    peaks <- unique(peaks)
    peaks$bdp_idx <- 1:length(peaks)
    anno <- annoPeaks(peaks, annoData, 
                      bindingType = "nearestBiDirectionalPromoters",
                      bindingRegion = c(-1*maxgap, maxgap))
    if(length(anno)<1){
        return(NA)
    }
    anno.s <- split(anno, anno$bdp_idx)
    len <- elementNROWS(anno.s)
    anno.s <- anno.s[len>=2]
    len <- sapply(anno.s, function(.ele){
        std <- .ele$feature.strand
        all(c("+", "-") %in% as.character(.ele$feature.strand))
    })
    anno.s <- anno.s[len]
    anno.s
}