test_that("getAllPeakSequence works not correct", {
    #Fixme, Zv9_scaffold.
    peaks <- GRanges(c("Zv9_scaffold3564","MT","1"),
                     IRanges(start=c(10,300,60348387),
                             end=c(50,330,60348389), 
                             names=c("p1","p2","p0")))
    suppressWarnings(seq <- 
                         getAllPeakSequence(peaks, 
                                            upstream = 20, 
                                            downstream = 20, 
                                            genome = Drerio))
    expect_equal(as.character(seq["p0",]$sequence),
                getSeq(Drerio, 
                       GRanges("chr1",
                               IRanges(60348367, 60348388)),
                       as.character=TRUE))
    expect_equal(as.character(seq["p1",]$sequence),
                suppressWarnings(getSeq(Drerio, 
                                        GRanges("Zv9_scaffold3564",
                                                IRanges(1, 70)),
                                        as.character=TRUE)))
    expect_equal(as.character(seq["p2",]$sequence),
                getSeq(Drerio, GRanges("chrM",
                                       IRanges(280, 350)),
                       as.character=TRUE))
})