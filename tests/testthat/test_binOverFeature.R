test_that("binOverFeature works not correct", {
    ##algorithm
    annotationData <- GRanges("X", IRanges(50, 50, name="X"), 
                              strand="+")
    peaks <- GRanges("X", IRanges(1:100, width=1), strand="+", 
                     score=rep(1, 100))
    bin <- binOverFeature(peaks, annotationData=annotationData, 
                          radius=50L, nbins=50L, minGeneLen=1L)
    expect_true(bin[1,1]==0)
    expect_equal(sum(bin[,1]), 99)
    
    peaks <- GRanges("X", IRanges(1:100, width=2), strand="+", 
                     score=rep(1, 100))
    bin <- binOverFeature(peaks, annotationData=annotationData, 
                          radius=50L, nbins=50L, minGeneLen=1L)
    expect_true(bin[1,1]==0)
    expect_true(bin[2,1]==1)
    expect_equal(sum(bin[,1]), 197)
    
    peaks <- GRanges("X", IRanges(1:50, 99:50), strand="+", score=rep(1, 50))
    bin <- binOverFeature(peaks, annotationData=annotationData, 
                          radius=50L, nbins=50L, minGeneLen=1L)
    expect_equal(sum(bin[,1]), sum(1:49)+sum(1:50))
    
    ## by simulation data
    annoGR <- annoGR(EnsDb.Hsapiens.v79)
    peaks <- as(annoGR, "GRanges")
    st <- as.character(strand(peaks))=="+"
    start(peaks)[st] <- start(peaks)[st]-100
    end(peaks)[st] <- start(peaks)[st]
    end(peaks)[!st] <- end(peaks)[!st]+100
    start(peaks)[!st] <- end(peaks)[!st]
    peaks$score <- 1
    bin <- binOverFeature(peaks, annotationData=annoGR, 
                          radius=500L, nbins=50L, minGeneLen=1L)
    expect_true(rownames(bin)[which(bin[,1]==max(bin[,1]))]=="-95")
})