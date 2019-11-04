test_that("toGRanges works not correct", {
    ##test format=c("BED", "GFF", "others")
    ##BED
    bed <- data.frame(cbind(space = c("1", "2"), start=c("100", "1000"),
                end=c("200", "1100"), name=c("peak1", "peak2")))
    gr <- toGRanges(bed, format="BED")
    expect_equal(start(gr), c(101, 1001))
    
    ##trim space
    bed$space <- paste(" ", 1:2)
    gr <- toGRanges(bed, format="BED")
    expect_equal(seqlevels(gr), c("1", "2"))
    bed$space <- paste(1:2, " ")
    gr <- toGRanges(bed, format="BED")
    expect_equal(seqlevels(gr), c("1", "2"))
    
    ##data.frame
    gr <- toGRanges(bed)
    expect_equal(start(gr), c(100, 1000))
    
    ##GFF
    GFF <- data.frame(cbind(space  = c("chr1", "chr2"), 
                            source=rep("Macs", 2),
                            feature=rep("peak", 2), 
                            start=c("100", "1000"), 
                            end=c("200", "1100"), 
                            score=c(60, 26),
                            strand=c(1, -1), 
                            frame=c(".", 2), 
                            group=c("peak1", "peak2")))
    gr <- toGRanges(GFF, format="GFF")
    
    ## EnsDb
    gr <- toGRanges(EnsDb.Hsapiens.v79, feature="gene")
    ## TxDb
    gr <- toGRanges(TxDb.Hsapiens.UCSC.hg38.knownGene, feature="gene")
})