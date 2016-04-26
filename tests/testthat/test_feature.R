test_that("featureAligend works not correct", {
    ##featureAlignedSignal
    cvglists <- list(A=RleList(chr1=Rle(rep(c(1:5, 5:1), 100))), 
                     B=RleList(chr1=Rle(rep(c(5:1, 1:5), 100))))
    feature.gr <- GRanges("chr1", IRanges(seq(15, 985, length.out=98), width=1))
    #feature.gr.shift <- shift(feature.gr, shift = -5)
    sig <- featureAlignedSignal(cvglists, feature.gr, 
                         upstream=5, downstream=5, 
                         n.tile=11)
    sig.unique <- lapply(sig, unique)
    expect_true(all(sapply(sig.unique, nrow)==1))
    expect_equal(sig.unique[["A"]][,1], 1)
    expect_equal(sig.unique[["B"]][,6], 1)
    f1 <- featureAlignedDistribution(sig, feature.gr, 
                               upstream=5, downstream=5, 
                               n.tile=11, type="l")
    f2 <- featureAlignedDistribution(cvglists, feature.gr,
                                     upstream=5, downstream=5, 
                                     n.tile=11, type="l")
    expect_identical(f1, f2)
    h1 <- featureAlignedHeatmap(sig, feature.gr,
                          upstream=5, downstream=5, 
                          n.tile=11, type="l")
    h2 <- featureAlignedHeatmap(cvglists, feature.gr,
                                upstream=5, downstream=5, 
                                n.tile=11, type="l")
    expect_identical(h1, h2)
    
    cvglists <- list(A=RleList(chr1=Rle(rep(c(0, 0, 3, 0, 0), 20))), 
                     B=RleList(chr1=Rle(rep(c(1:5, 5:1), 10))))
    feature.gr <- GRanges("chr1", IRanges(seq(15, 85, length.out=8), width=1))
    sig <- featureAlignedSignal(cvglists, feature.gr, 
                                upstream=4, downstream=5, 
                                n.tile=10)
    f<- featureAlignedDistribution(cvglists, feature.gr,
                               upstream=4, downstream=5, 
                               n.tile=10, type="l")
    expect_equal(f[, "A"], c(0, 0, 3, 0, 0, 0, 0, 3, 0, 0))
    expect_equal(f[, "B"], c(1:5, 5:1))
})