test_annoPeaks <- function(){
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
    feature <- annoGR(feature)
    anno <- lapply(list(c(-1,1), c(-2,1), c(-2,3), c(-3,3), c(-4,4)), 
                   annoPeaks, peaks=myPeak, annoData=feature, 
                   bindingType="startSite", ignore.peak.strand=TRUE)
    anno.name <- sapply(anno, function(.ele) 
        paste0(names(.ele), .ele$feature))
    checkEquals(anno.name[[1]], "aA")
    checkEquals(anno.name[[2]], c("aA", "aD", "aG"))
    checkEquals(anno.name[[3]], c("aA", "aD", "aG"))
    checkEquals(anno.name[[4]], c("aA", "aD", "aF", "aG", "aK"))
    checkEquals(anno.name[[5]], 
                c("aA", "aD", "aF", "aG", 'aJ', "aK", "aM", "aN"))
    anno <- lapply(c("startSite", "endSite", "fullRange", "bothSidesNearest"),
                    annoPeaks, peaks=myPeak, annoData=feature, 
                    bindingRegion=c(-3, 3), ignore.peak.strand=TRUE)
    anno.name <- sapply(anno, function(.ele) 
        paste0(names(.ele), .ele$feature))
    checkEquals(anno.name[[1]], c("aA", "aD", "aF", "aG", "aK"))
    checkEquals(anno.name[[2]], c("aB", "aC", "aE", "aM", "aN"))
    checkEquals(anno.name[[3]], c("aA", "aB", "aC", "aD", "aE", 
                                  "aF", "aG", "aK", "aM", "aN"))
    checkEquals(anno.name[[4]], c("aA", "aB", "aC", "aD", "aM", "aN"))
}