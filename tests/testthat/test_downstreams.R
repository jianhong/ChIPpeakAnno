test_that("downstreams works not correct", {
    ##GRanges
    ## myPeak
    ## 1 2 3 4 5 6 7 8 9 10
    ##       a a a b b       + -
    
    myPeak = GRanges("1", 
                     IRanges(start=c(4, 7), 
                             end=c(6, 8), 
                             names=letters[1:2]), 
                     strand=c("+", "-"))
    dws <- downstreams(myPeak, upstream = 0, downstream = 2)
    expect_equal(width(dws), c(2, 2))
    expect_equal(start(dws), c(7, 5))
    expect_equal(strand(dws), strand(myPeak))
    ##0 bdp
    dws <- downstreams(myPeak, upstream = 0, downstream = 0)
    expect_equal(width(dws), c(0, 0))
    ##up
    dws <- downstreams(myPeak, upstream = 2, downstream = 0)
    expect_equal(width(dws), c(2, 2))
    expect_equal(start(dws), c(5, 7))
    expect_equal(end(dws), end(myPeak))
    expect_equal(strand(dws), strand(myPeak))
    ##up & down
    dws <- downstreams(myPeak, upstream = 2, downstream = 2)
    expect_equal(width(dws), c(4, 4))
    expect_equal(start(dws), c(5, 5))
    expect_equal(end(dws), c(8, 8))
    expect_equal(strand(dws), strand(myPeak))
})