test_that("annoPeaks works not correct", {
    #check start
    #               a
    #      12345678901234567890
    #            -AA -BB
    #            +CC +DD
    #           +EE   +FF
    #           -GG   -HH
    #          +II     +JJ
    #          -KK     -LL
    #           +MMMMMMM
    #           -NNNNNNN
    myPeak = GRanges("chr1", 
                     IRanges(start=10, 
                             end=10, 
                             names="a"), 
                     strand="*")
    feature = GRanges("chr1", 
                      IRanges(start=c(8, 12, 8, 12, 7, 13, 7, 
                                      13, 6, 14, 6, 14, 7, 7),
                              end=c(9, 13, 9, 13, 8, 14, 8, 
                                    14, 7, 15, 7, 15, 13, 13),
                              names=LETTERS[1:14]), 
                      strand=c("-", "-", "+", "+", "+", "+", "-", 
                               "-", "+", "+", "-", "-", "+", "-"))
    expect_error(annoPeaks(myPeak, feature, bindingRegion=c(10, 10)))
    feature <- annoGR(feature)
    anno <- lapply(list(c(-1,1), c(-2,1), c(-2,3), c(-3,3), c(-4,4)), 
                   annoPeaks, peaks=myPeak, annoData=feature, 
                   bindingType="startSite", ignore.peak.strand=TRUE)
    anno.name <- sapply(anno, function(.ele) 
        paste0(names(.ele), .ele$feature))
    expect_equal(anno.name[[1]], "aA")
    expect_equal(anno.name[[2]], c("aA", "aD", "aG"))
    expect_equal(anno.name[[3]], c("aA", "aD", "aG"))
    expect_equal(anno.name[[4]], c("aA", "aD", "aF", "aG", "aK"))
    expect_equal(anno.name[[5]], 
                c("aA", "aD", "aF", "aG", 'aJ', "aK", "aM", "aN"))
    anno <- lapply(c("startSite", "endSite", "fullRange", "bothSidesNearest"),
                    annoPeaks, peaks=myPeak, annoData=feature, 
                    bindingRegion=c(-3, 3), ignore.peak.strand=TRUE)
    anno.name <- sapply(anno, function(.ele) 
        paste0(names(.ele), .ele$feature))
    expect_equal(anno.name[[1]], c("aA", "aD", "aF", "aG", "aK"))
    expect_equal(anno.name[[2]], c("aB", "aC", "aE", "aM", "aN"))
    expect_equal(anno.name[[3]], c("aA", "aB", "aC", "aD", "aE", 
                                  "aF", "aG", "aK", "aM", "aN"))
    expect_equal(anno.name[[4]], c("aA", "aB", "aC", "aD", "aM", "aN"))
    
    #             AA
    #      12345678901234567890
    #     -aaaa      +cccc
    #       +bbbb
    annoData<- GRanges(seqnames = "chr1", 
                       ranges=IRanges(start=c(1,3,12), 
                                      end=c(4,6,15),
                                      names=letters[1:3]), 
                       strand = c("-", "+", "+"))
    peak<- GRanges(seqnames = "chr1", 
                   ranges=IRanges(start=8, width=2), 
                   strand = "+")
    peak.anno <- annoPeaks(peak, annoData)
    expect_equal(peak.anno$feature, c("a", "c"))
    peak.anno <- annoPeaks(peak, annoData, 
                           bindingType="nearestBiDirectionalPromoters",
                           bindingRegion=c(-2, 1), 
                           ignore.peak.strand=TRUE)
    expect_true(length(peak.anno)==0)
    peak.anno <- annoPeaks(peak, annoData, 
                           bindingType="nearestBiDirectionalPromoters",
                           bindingRegion=c(-3, 1), 
                           ignore.peak.strand=TRUE)
    expect_equal(peak.anno$feature, "c")
    #             AA
    #      12345678901234567890
    #     -aaaa      +cccc
    #       +bbbbbb
    #       +ddd
    annoData<- GRanges(seqnames = "chr1", 
                       ranges=IRanges(start=c(1,3,12,3), 
                                      end=c(4,8,15,5),
                                      names=letters[1:4]), 
                       strand = c("-", "+", "+", "+"))
    peak.anno <- annoPeaks(peak, annoData, 
                           bindingType="nearestBiDirectionalPromoters",
                           bindingRegion=c(-3, 1), 
                           ignore.peak.strand=TRUE)
    expect_equal(peak.anno$feature, "c")
    peak.anno <- annoPeaks(peak, annoData, 
                           bindingType="nearestBiDirectionalPromoters",
                           bindingRegion=c(-4, 1), 
                           ignore.peak.strand=TRUE)
    expect_equal(peak.anno$feature, c("a", "c"))
    peak.anno <- annoPeaks(peak, annoData, 
                           bindingType="nearestBiDirectionalPromoters",
                           bindingRegion=c(-4, 6), 
                           ignore.peak.strand=TRUE)
    expect_equal(peak.anno$feature, c("a", "c"))
    peak.anno <- annoPeaks(peak, annoData, 
                           bindingType="nearestBiDirectionalPromoters",
                           bindingRegion=c(-3, 6), 
                           ignore.peak.strand=TRUE)
    expect_equal(peak.anno$feature, "c")
})