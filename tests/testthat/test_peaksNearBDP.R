test_that("peaksNearBDP works not correct", {
    ##GRanges
    ## myPeak
    ## 1 2 3 4 5 6 7 8 9 10
    ##   a a a a a b b   c
    ## A A B         D D    +
    ##     F F F   C C      -
    ##   E E E E E          +
    ## G G G G G            +
    ##     H H H H H        +
    ## I I I I I I I        -
    
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
    expect_equal(annotatedBDP$n.peaksWithBDP, 2)
    ##0 bdp
    strand(feature) <- "+"
    suppressWarnings(annotatedBDP <- 
                         peaksNearBDP(myPeak, 
                                      AnnotationData=feature,
                                      MaxDistance=2,
                                      PeakLocForDistance =  "middle", 
                                      FeatureLocForDistance = "TSS"))
    expect_equal(annotatedBDP$n.peaksWithBDP, 0)
})