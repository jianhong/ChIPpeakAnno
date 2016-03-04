test_that("reCenterPeaks works not correct", {
    ## zero length
    reCenterPeaks(GRanges(), 1)
})