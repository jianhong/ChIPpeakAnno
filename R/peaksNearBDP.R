peaksNearBDP <- function(myPeakList, mart, AnnotationData, MaxDistance=5000,PeakLocForDistance = c("start", "middle", "end"), 
FeatureLocForDistance = c("TSS", "middle","start", "end","geneEnd"))
{
	if (missing(myPeakList)) {
        stop("Missing required argument myPeakList!")
    }
    if (class(myPeakList) != "RangedData") {
        stop("No valid myPeakList passed in. It needs to be RangedData object")
    }
    PeakLocForDistance = match.arg(PeakLocForDistance)
    FeatureLocForDistance = match.arg(FeatureLocForDistance)
    if (!missing(AnnotationData))
	{
		AnnoPlus =AnnotationData[AnnotationData$strand =="+"  | AnnotationData$strand ==1,]
		AnnoMinus =AnnotationData[AnnotationData$strand =="-" | AnnotationData$strand ==-1,]
	}
	else if (missing(mart) || class(mart) != "Mart") {
          		  stop("Error in querying biomart database. No valid mart object is passed in! Suggest call getAnnotation before calling PeaksNearBDP")
      	  }
	else
	{
       	AnnotationData <- getAnnotation(mart, feature = "TSS") 
		message("Done querying biomart database, start annotating ....Better way would be calling getAnnotation before calling peaksNearBDP")
		AnnoPlus =AnnotationData[AnnotationData$strand =="+"  | AnnotationData$strand ==1,]
		AnnoMinus =AnnotationData[AnnotationData$strand =="-" | AnnotationData$strand ==-1,]
	}
    plus = annotatePeakInBatch(myPeakList, AnnotationData = AnnoPlus,PeakLocForDistance=PeakLocForDistance, FeatureLocForDistance=FeatureLocForDistance)
	minus = annotatePeakInBatch(myPeakList, AnnotationData = AnnoMinus)
	plus.passed = plus[!is.na(plus$shortestDistance) & plus$shortestDistance<=MaxDistance, ]
	minus.passed = minus[!is.na(minus$shortestDistance) & minus$shortestDistance <=MaxDistance, ]
	passed = rbind(plus.passed, minus.passed)
	summary = as.data.frame(table(passed$peak))
	BDP = summary[summary[,2] >1,]
	peaksWithBDP = passed[passed$peak %in% BDP[,1],]
	temp = as.data.frame(table(summary[,2]))
	percentPeaksWithBDP  = temp[2,2]/dim(myPeakList)[1]
	list(peaksWithBDP = peaksWithBDP, percentPeaksWithBDP = percentPeaksWithBDP, n.peaks=dim(myPeakList)[1], n.peaksWithBDP=temp[2,2])
}
