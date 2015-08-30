test_peaksNearBDP<-function(){
    ##GRanges
    myPeak = GRanges("1", 
                     IRanges(start=c(2, 7, 10), 
                             end=c(6, 8, 10), 
                             names=letters[1:3]), 
                     strand="*")
    feature = GRanges("1", 
                      IRanges(start=c(1,3, 6, 8, 2, 3, 1, 3, 1), 
                              end=c(2, 3, 7, 9, 6, 5, 5, 7, 7), 
                              names=LETTERS[1:9]), 
                      strand=c("+","+", "-", "+", "+", "-", "+", "+", "-"))
    annotatedBDP = peaksNearBDP(myPeak, AnnotationData=feature,
                                MaxDistance=2,PeakLocForDistance =  "middle", 
                                FeatureLocForDistance = "TSS")
    checkEquals(annotatedBDP$n.peaksWithBDP, 2)
    ##0 bdp
    strand(feature) <- "+"
    suppressWarnings(annotatedBDP <- 
                         peaksNearBDP(myPeak, 
                                      AnnotationData=feature,
                                      MaxDistance=2,
                                      PeakLocForDistance =  "middle", 
                                      FeatureLocForDistance = "TSS"))
    checkEquals(annotatedBDP$n.peaksWithBDP, 0)
}