test_that("findEnhancers works not correct", {
    ## expect GRanges()
    DNAinteractiveData <- GRanges("chr1", IRanges(100, 1500), 
                                  blocks=IRangesList(IRanges(c(1, 901), 
                                                             c(101, 1401))))
    peaks <- GRanges("chr1", IRanges(c(50, 150, 250, 950, 1050, 1450, 1550), 
                                    width=10, names=letters[1:7]))
    anno <- GRanges("chr1", IRanges(c(40, 1560), width=10, names=LETTERS[1:2]),
                    strand=c("-", "+"))
    enhancers <- findEnhancers(peaks, anno, DNAinteractiveData,
                                bindingType="nearestBiDirectionalPromoters",
                                bindingRegion=c(-120, 120),
                                ignore.peak.strand = TRUE)
    expect_equal(letters[1:7], enhancers$peak)
    expect_equal(c(100, 100, 109, 91, 100, 91, 100), enhancers$distance)
})