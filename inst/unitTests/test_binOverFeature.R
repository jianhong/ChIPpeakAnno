test_binOverFeature<-function(){
    ##algorithm
    annotationData <- GRanges("x", IRanges(50, 50, name="X"), 
                              strand="+")
    peaks <- GRanges("x", IRanges(1:100, width=1), strand="+", 
                     score=rep(1, 100))
    bin <- binOverFeature(peaks, annotationData=annotationData, 
                          radius=50L, nbins=50L, minGeneLen=1L)
    checkTrue(bin[1,1]==0)
    checkEqualsNumeric(sum(bin[,1]), 99)
    
    peaks <- GRanges("x", IRanges(1:100, width=2), strand="+", 
                     score=rep(1, 100))
    bin <- binOverFeature(peaks, annotationData=annotationData, 
                          radius=50L, nbins=50L, minGeneLen=1L)
    checkTrue(bin[1,1]==0)
    checkTrue(bin[2,1]==1)
    checkEqualsNumeric(sum(bin[,1]), 197)
    
    peaks <- GRanges("x", IRanges(1:50, 99:50), strand="+", score=rep(1, 50))
    bin <- binOverFeature(peaks, annotationData=annotationData, 
                          radius=50L, nbins=50L, minGeneLen=1L)
    checkEqualsNumeric(sum(bin[,1]), sum(1:49)+sum(1:50))
}