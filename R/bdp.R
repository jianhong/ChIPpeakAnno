#' obtain the peaks near bi-directional promoters
#' 
#' Obtain the peaks near bi-directional promoters. Also output percent of peaks
#' near bi-directional promoters.
#' 
#' 
#' @param peaks peak list, \link[GenomicRanges:GRanges-class]{GRanges} object
#' @param annoData annotation data, \link{annoGR} object
#' @param maxgap maxgap between peak and TSS
#' @param ...  Not used.
#' @return Output is a list of GRanges object of the peaks near bi-directional
#' promoters.
#' @author Jianhong Ou
#' @seealso See Also as \code{\link{annoPeaks}}, \code{\link{annoGR}}
#' @keywords misc
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom S4Vectors elementNROWS
#' @examples
#' 
#'   if(interactive() || Sys.getenv("USER")=="jianhongou"){
#'     library(ensembldb)
#'     library(EnsDb.Hsapiens.v75)
#'     data("myPeakList")
#'     annoGR <- annoGR(EnsDb.Hsapiens.v75)
#'     seqlevelsStyle(myPeakList) <- seqlevelsStyle(annoGR)
#'     ChIPpeakAnno:::bdp(myPeakList, annoGR)
#'   }
#' 
bdp <- function(peaks, annoData, maxgap=2000L, ...){
    stopifnot(inherits(peaks, "GRanges"))
    stopifnot(inherits(annoData, c("annoGR", "GRanges")))
    stopifnot(length(intersect(seqlevelsStyle(peaks),seqlevelsStyle(annoData)))>0)
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
