test_featureAligend<-function(){
    ##featureAlignedSignal
    cvglists <- list(A=RleList(chr1=Rle(rep(c(1:5, 5:1), 100))), 
                     B=RleList(chr1=Rle(rep(c(5:1, 1:5), 100))))
    feature.gr <- GRanges("chr1", IRanges(seq(15, 985, length.out=98), width=1))
    sig <- featureAlignedSignal(cvglists, feature.gr, 
                         upstream=5, downstream=5, 
                         n.tile=11)
    sig.unique <- lapply(sig, unique)
    checkTrue(all(sapply(sig.unique, nrow)==1))
    checkEquals(sig.unique[["A"]][,1], 1)
    checkEquals(sig.unique[["B"]][,6], 1)
    f1 <- featureAlignedDistribution(sig, feature.gr, 
                               upstream=5, downstream=5, 
                               n.tile=11, type="l")
    f2 <- featureAlignedDistribution(cvglists, feature.gr,
                                     upstream=5, downstream=5, 
                                     n.tile=11, type="l")
    checkIdentical(f1, f2)
    h1 <- featureAlignedHeatmap(sig, feature.gr,
                          upstream=5, downstream=5, 
                          n.tile=11, type="l")
    h2 <- featureAlignedHeatmap(cvglists, feature.gr,
                                upstream=5, downstream=5, 
                                n.tile=11, type="l")
    checkIdentical(h1, h2)
}