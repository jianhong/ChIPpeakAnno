test_that("toGRanges works not correct", {
    ##test format=c("BED", "GFF", "RangedData", "others")
    ##RangedData
    if(exists("RangedData")){
        rd <- RangedData(IRanges(start = c(967654, 2010897, 2496704),
                                 end = c(967754, 2010997, 2496804), 
                                 names = c("Site1", "Site2", "Site3")),
                         space = c("1", "2", "3"), strand=as.integer(1),
                         feature=c("a","b","f"))
        gr <- toGRanges(rd, format="RangedData")
        expect_equal(start(rd), start(gr))
        expect_equal(end(rd), end(gr))
    }
    
    ##BED
    bed <- data.frame(cbind(space = c("1", "2"), start=c("100", "1000"),
                end=c("200", "1100"), name=c("peak1", "peak2")))
    gr <- toGRanges(bed, format="BED")
    expect_equal(start(gr), c(101, 1001))
    
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
})