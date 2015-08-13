buildBindingDistribution <- function(x, AnnotationData,
                                     bindingType=c("TSS", "geneEnd"), 
                                     featureType=c("transcript", "exon")){
    if(class(x)!="GRanges"){
        stop("x must be an object of GRanges")
    }
    if(class(AnnotationData)!="GRanges"){
        stop("AnnotationData must be an object of GRanges")
    }
    bindingType <- match.arg(bindingType)
    featureType <- match.arg(featureType)
    AnnotationData <- unique(AnnotationData)
    suppressWarnings({
        anno <- annotatePeakInBatch(x, AnnotationData=AnnotationData, 
                                    output="shortestDistance",
                                    FeatureLocForDistance=bindingType, 
                                    PeakLocForDistance="middle")
    })
    
    anno <- unique(anno)
    if(length(anno)<1){
        stop("Can not annotate input x with the AnnotationData.")
    }
    distancetoFeature <- anno$distancetoFeature
    h <- hist(distancetoFeature, breaks=ceiling(length(x)/100), 
              plot=FALSE)
    N <- as.integer(h$counts)
    diff <- as.integer(floor(diff(h$breaks)[1]/2))
    offset <- as.integer(floor(h$mids))
    offset <- offset[N>0]
    N <- N[N>0]
    new("bindist", counts=N,
        mids=offset,
        halfBinSize=diff,
        bindingType=bindingType,
        featureType=featureType)
}