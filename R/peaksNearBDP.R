#' obtain the peaks near bi-directional promoters
#' 
#' Obtain the peaks near bi-directional promoters. Also output percent of peaks
#' near bi-directional promoters.
#' 
#' 
#' @param myPeakList \link[GenomicRanges:GRanges-class]{GRanges}: See example
#' below
#' @param AnnotationData annotation data obtained from getAnnotation or
#' customized annotation of class \link[GenomicRanges:GRanges-class]{GRanges}
#' containing additional variable: strand (1 or + for plus strand and -1 or -
#' for minus strand). For example, data(TSS.human.NCBI36),
#' data(TSS.mouse.NCBIM37), data(TSS.rat.RGSC3.4) and data(TSS.zebrafish.Zv8).
#' @param MaxDistance Specify the maximum gap allowed between the peak and
#' nearest gene
#' @param ... Not used
#' @return A list of 4 \item{list("peaksWithBDP")}{annotated Peaks containing
#' bi-directional promoters.
#' 
#' GRangesList with slot start holding the start position of the peak, slot end
#' holding the end position of the peak, slot space holding the chromosome
#' location where the peak is located, slot rownames holding the id of the
#' peak.  In addition, the following variables are included.
#' 
#' feature: id of the feature such as ensembl gene ID
#' 
#' insideFeature: upstream: peak resides upstream of the feature; downstream:
#' peak resides downstream of the feature; inside: peak resides inside the
#' feature; overlapStart: peak overlaps with the start of the feature;
#' overlapEnd: peak overlaps with the end of the feature; includeFeature: peak
#' include the feature entirely.
#' 
#' distancetoFeature: distance to the nearest feature such as transcription
#' start site.  By default, the distance is calculated as the distance between
#' the start of the binding site and the TSS that is the gene start for genes
#' located on the forward strand and the gene end for genes located on the
#' reverse strand. The user can specify the location of peak and location of
#' feature for calculating this
#' 
#' feature_range: start and end position of the feature such as gene
#' 
#' feature_strand: 1 or + for positive strand and -1 or - for negative strand
#' where the feature is located } \item{list("percentPeaksWithBDP")}{The
#' percent of input peaks containing bi-directional promoters}
#' \item{list("n.peaks")}{The total number of input peaks}
#' \item{list("n.peaksWithBDP")}{The # of input peaks containing bi-directional
#' promoters}
#' @author Lihua Julie Zhu, Jianhong Ou
#' @seealso annotatePeakInBatch, findOverlappingPeaks, makeVennDiagram
#' @references Zhu L.J. et al. (2010) ChIPpeakAnno: a Bioconductor package to
#' annotate ChIP-seq and ChIP-chip data. BMC Bioinformatics 2010,
#' 11:237doi:10.1186/1471-2105-11-237
#' @keywords misc
#' @export
#' @examples
#' 
#' 
#' if (interactive() || Sys.getenv("USER")=="jianhongou")
#' {
#'     data(myPeakList)
#'     data(TSS.human.NCBI36)
#'     seqlevelsStyle(TSS.human.NCBI36) <- seqlevelsStyle(myPeakList)
#'     annotatedBDP = peaksNearBDP(myPeakList[1:6,], 
#'                                 AnnotationData=TSS.human.NCBI36,
#'                                 MaxDistance=5000,
#'                                 PeakLocForDistance =  "middle", 
#'                                 FeatureLocForDistance = "TSS")
#'     c(annotatedBDP$percentPeaksWithBDP, annotatedBDP$n.peaks, 
#'       annotatedBDP$n.peaksWithBDP)
#' }
#' 
peaksNearBDP <- function(myPeakList, AnnotationData,
                         MaxDistance=5000L, ...){
    if (missing(myPeakList)) {
        stop("Missing required argument myPeakList!")
    }
    if (!inherits(myPeakList, c("GRanges"))) {
        stop("myPeakList needs to be GRanges object")
    }
    if (!missing(AnnotationData)){        
        if (!inherits(AnnotationData, c("GRanges", "annoGR"))) {
            stop("AnnotationData needs to be GRanges or annoGR object")
        }
        if(is(AnnotationData, "annoGR"))
            AnnotationData <- AnnotationData@gr
    }else{
        stop("Missing required argument AnnotationData!")
    }
    stopifnot(length(intersect(seqlevelsStyle(myPeakList),
                               seqlevelsStyle(AnnotationData)))>0)
    stopifnot(is.numeric(MaxDistance))
    
    MaxDistance <- round(MaxDistance[1])
    myPeakList <- unique(myPeakList)
    myPeakList$bdp_idx <- seq_along(myPeakList)
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
