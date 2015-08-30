peaksNearBDP <- function(myPeakList, mart, AnnotationData, 
                         MaxDistance=5000,
                         PeakLocForDistance = c("start", "middle", "end"), 
                         FeatureLocForDistance = c("TSS", "middle","start", 
                                                   "end","geneEnd")){
    if (missing(myPeakList)) {
        stop("Missing required argument myPeakList!")
    }
    if (!inherits(myPeakList, c("RangedData", "GRanges"))) {
        stop("myPeakList needs to be GRanges object")
    }
    PeakLocForDistance = match.arg(PeakLocForDistance)
    FeatureLocForDistance = match.arg(FeatureLocForDistance)
    if (!missing(AnnotationData))
    {        
        if (!inherits(AnnotationData, c("RangedData", "GRanges", "annoGR"))) {
            stop("AnnotationData needs to be GRanges or annoGR object")
        }
        if(class(AnnotationData)=="annoGR")
            AnnotationData <- AnnotationData@gr
        if(class(AnnotationData)=="RangedData") 
            AnnotationData <- toGRanges(AnnotationData, format="RangedData")
    }
    else if (missing(mart) || class(mart) != "Mart") {
        stop("Error in querying biomart database. 
             No valid mart object is passed in! 
             Suggest call getAnnotation before calling PeaksNearBDP")
    }
    else
    {
        AnnotationData <- getAnnotation(mart, featureType = "TSS")
        message("Done querying biomart database, start annotating ....
                Better way would be calling getAnnotation before 
                calling peaksNearBDP")
    }
    AnnoPlus = AnnotationData[strand(AnnotationData) =="+" | 
                                  strand(AnnotationData) ==1]
    AnnoMinus =AnnotationData[strand(AnnotationData) =="-" | 
                                  strand(AnnotationData) ==-1]
    if(class(myPeakList)=="RangedData") 
        myPeakList <- toGRanges(myPeakList, format="RangedData")
    plus = annotatePeakInBatch(myPeakList, AnnotationData = AnnoPlus,
                               PeakLocForDistance=PeakLocForDistance, 
                               FeatureLocForDistance=FeatureLocForDistance)
    minus = annotatePeakInBatch(myPeakList, AnnotationData = AnnoMinus, 
                                PeakLocForDistance=PeakLocForDistance, 
                                FeatureLocForDistance=FeatureLocForDistance)
    plus.passed = plus[!is.na(plus$shortestDistance) & 
                           plus$shortestDistance<=MaxDistance]
    minus.passed = minus[!is.na(minus$shortestDistance) & 
                             minus$shortestDistance <=MaxDistance]
    passed = c(plus.passed, minus.passed)
    summary = 
        table(as.data.frame(mcols(passed)[, c("peak", "feature_strand")]))
    if(nrow(summary)<1) return(NULL)
    isBDP <- apply(summary, 1, function(.ele) sum(.ele>0)) > 1
    BDP = summary[isBDP, , drop=FALSE]
    peaksWithBDP = passed[passed$peak %in% rownames(BDP)]
    peaksWithBDP <- peaksWithBDP[order(as.character(seqnames(peaksWithBDP)), 
                                       peaksWithBDP$peak, 
                                       peaksWithBDP$feature_strand)]
    temp = nrow(BDP)
    percentPeaksWithBDP  = temp/length(myPeakList)
    list(peaksWithBDP = peaksWithBDP, 
         percentPeaksWithBDP = percentPeaksWithBDP, 
         n.peaks=length(myPeakList), 
         n.peaksWithBDP=temp)
}
