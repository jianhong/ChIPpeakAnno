test_that("filterByOverlaps works not correct", {
    ##GRanges
    ## myPeak
    ## 1 2 3 4 5 6 7 8 9 10
    ##       a a a a         +
    ##           b b b       -
    
    myPeak = GRanges("1", 
                     IRanges(start=c(4, 6), 
                             end=c(7, 8), 
                             names=letters[1:2]), 
                     strand=c("+", "-"))
    o <- filterByOverlaps(myPeak[1], myPeak[2], ignore.strand=FALSE)
    expect_identical(o, unname(myPeak[1]))
    o <- filterByOverlaps(myPeak[1], myPeak[2], ignore.strand=TRUE)
    expect_equal(end(o), 5)
})
