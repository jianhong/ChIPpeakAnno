#' @importFrom graphics hist
#' @import methods
#' @importFrom grDevices boxplot.stats
buildBindingDistribution <- function(x, AnnotationData,
                                     bindingType=c("TSS", "geneEnd"), 
                                     featureType=c("transcript", "exon")){
    if(missing(AnnotationData) || missing(x)){
        stop("x and AnnotationData is required.")
    }
    if(!is(x, "GRanges")){
        stop("x must be an object of GRanges")
    }
    if(class(AnnotationData)=="annoGR"){
        AnnotationData <- as(AnnotationData, "GRanges")
    }
    if(!is(AnnotationData, "GRanges")){
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
    ## calculate the breaks
    ## previous just use ceiling(length(x)/100)
    ## now use boxplot.stats to remove the out
    box.stats <- boxplot.stats(distancetoFeature, coef=1.5, 
                               do.conf=FALSE, do.out=FALSE)
    whiskers <- box.stats$stats[5] - box.stats$stats[1]
    minWid <- 100
    if(whiskers>100*minWid){
        breaks <- seq(box.stats$stats[1], box.stats$stats[5],
                      length.out=100)
        minWid <- diff(breaks)[1]
    }else{
        if(whiskers>minWid){
            breaks <- seq(floor(box.stats$stats[1]/minWid), 
                          ceiling(box.stats$stats[5]/minWid),
                          by=1)*minWid
        }else{
            breaks <- mean(box.stats$stats[c(1,5)])
            breaks <- c(breaks-minWid/2, breaks+minWid/2)
        }
    }
    if(min(distancetoFeature)<min(breaks)) 
        breaks <- c(min(distancetoFeature), breaks)
    if(max(distancetoFeature)>max(breaks))
        breaks <- c(breaks, max(distancetoFeature))
    h <- hist(distancetoFeature, breaks=breaks, 
              plot=FALSE)
    N <- as.integer(h$counts)
    diff <- as.integer(floor(minWid/2))
    offset <- as.integer(floor(h$mids))
    offset <- offset[N>0]
    N <- N[N>0]
    new("bindist", counts=N,
        mids=offset,
        halfBinSize=diff,
        bindingType=bindingType,
        featureType=featureType)
}