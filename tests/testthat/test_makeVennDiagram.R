test_that("makeVennDiagram works not correct", {
    peaks1 = GRanges(seqnames = c("1", "2", "3"), 
                     IRanges(start = c(967654, 2010897, 2496704),
                             end = c(967754, 2010997, 2496804), 
                             names = c("Site1", "Site2", "Site3")),
                     strand="+",
                     feature=c("a","b","f"))
    peaks2 = GRanges(seqnames = c("1", "2", "3", "1", "2"), 
                     IRanges(start = c(967659, 2010898, 2496700,
                                       3075866, 3123260),
                             end = c(967869, 2011108, 2496920, 
                                     3076166, 3123470),
                             names = c("t1", "t2", "t3", "t4", "t5")), 
                     strand = c("+", "+", "+","-", "+"), 
                     feature=c("a","b","c","d","a"))
    x<-makeVennDiagram(GRangesList(peaks1, peaks2), 
                       NameOfPeaks=c("TF1", "TF2"),
                       totalTest=100)
    expect_true(x$p.value[1,"pval"] - 6.184292e-05 < 1e-05)
    expect_equal(x$vennCounts[,"Counts"], c(95, 2, 0, 3))
    
    y<-makeVennDiagram(GRangesList(peaks1, peaks2, peaks1, peaks2), 
                       NameOfPeaks=c("TF1", "TF2","TF3", "TF4"),
                       totalTest=100, by="feature")
    expect_true(y$p.value[y$p.value[, "TF1"]==1 & 
                                     y$p.value[, "TF2"]==1,"pval"] - 
                       0.003586889 < 0.001)
    expect_true(y$p.value[y$p.value[, "TF1"]==1 & 
                                     y$p.value[, "TF3"]==1,"pval"] - 
                       6.184292e-06 < 1e-06)
    expect_true(y$p.value[y$p.value[, "TF2"]==1 & 
                                     y$p.value[, "TF3"]==1,"pval"] -
                       0.003586889 < 0.001)
    expect_true(y$p.value[y$p.value[, "TF1"]==1 & 
                                     y$p.value[, "TF4"]==1,"pval"] - 
                       0.003586889 < 0.001)
    expect_true(y$p.value[y$p.value[, "TF3"]==1 & 
                                     y$p.value[, "TF4"]==1,"pval"] -
                       0.003586889 < 0.001)
    expect_equal(y$vennCounts[,"Counts"], 
                       c(95, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2))
    ##test 0 overlaps
    peaks1 = GRanges(seqnames = c("1", "1", "1"), 
                     IRanges(start = c(1, 100, 200),
                             end = c(2, 101, 201), 
                             names = c("Site1", "Site2", "Site3")),
                     strand="+",
                     feature=c("a","b","f"))
    peaks2 = GRanges(seqnames = c("1", "1", "1"), 
                     IRanges(start = c(1, 100, 200),
                             end = c(2, 101, 201), 
                             names = c("Site1", "Site2", "Site3")),
                     strand="+",
                     feature=c("a","b","f"))
    suppressWarnings(makeVennDiagram(GRangesList(peaks1,peaks2), 
                                     NameOfPeaks=c("TF1", "TF2"),
                                     totalTest=100,scaled=F, euler.d=F))
    
    peaks1 = GRanges(seqnames = c("1", "1", "1"), 
                     IRanges(start = c(1, 100, 200),
                             end = c(2, 101, 201), 
                             names = c("Site1", "Site2", "Site3")),
                     strand="+",
                     feature=c("a","b","f"))
    peaks2 = GRanges(seqnames = c("2", "2", "2"), 
                     IRanges(start = c(1, 100, 200),
                             end = c(2, 101, 201), 
                             names = c("Site1", "Site2", "Site3")),
                     strand="+",
                     feature=c("a","b","f"))
    suppressWarnings(makeVennDiagram(GRangesList(peaks1,peaks2), 
                                     NameOfPeaks=c("TF1", "TF2"),
                                     totalTest=100,scaled=F, euler.d=F))
    
    peaks1 = GRanges(seqnames = c("1", "1", "1"), 
                     IRanges(start = c(1, 100, 200),
                             end = c(2, 101, 201), 
                             names = c("Site1", "Site2", "Site3")),
                     strand="+",
                     feature=c("a","b","f"))
    peaks2 = GRanges(seqnames = c("1", "1", "1"),
                     IRanges(start = c(5, 105, 205),
                             end = c(6, 106, 206), 
                             names = c("Site1", "Site2", "Site3")), 
                     strand="+",
                     feature=c("a","b","f"))
    suppressWarnings(makeVennDiagram(GRangesList(peaks1,peaks2), 
                                     NameOfPeaks=c("TF1", "TF2"),
                                     totalTest=100,scaled=F, euler.d=F))
    
    peaks1 = GRanges(seqnames = c("1", "1", "1"),
                     IRanges(start = c(1, 100, 200),
                             end = c(2, 101, 201), 
                             names = c("Site1", "Site2", "Site3")),
                     strand="+",
                     feature=c("a","b","f"))
    peaks2 = GRanges(IRanges(start = c(1, 100, 200),
                             end = c(1, 101, 201), 
                             names = c("Site1", "Site2", "Site3")),
                     seqnames = c("2", "2", "2"),
                     strand="+",feature=c("a","b","f"))
    peaks3 = GRanges(seqnames = c("1", "2", "2"),
                     IRanges(start = c(1, 100, 200),
                             end = c(3, 103, 203), 
                             names = c("Site1", "Site2", "Site3")),
                     strand="+",feature=c("a","b","f"))
    peaks4 = GRanges(seqnames = c("1", "2", "2"),
                     IRanges(start = c(1, 100, 200),
                             end = c(2, 101, 201), 
                             names = c("Site1", "Site2", "Site3")),
                     strand="+",
                     feature=c("a","b","f"))
    suppressWarnings(makeVennDiagram(GRangesList(peaks1,peaks2,peaks3),
                                     NameOfPeaks=c("TF1", "TF2","TF3"),
                                     totalTest=100,scaled=F, euler.d=F))
    #data same as above except first chrom of first range in peaks2
    peaks1 = GRanges(seqnames = c("1", "1", "1"),
                     IRanges(start = c(1, 100, 200),
                             end = c(2, 101, 201), 
                             names = c("Site1", "Site2", "Site3")),
                     strand="+",feature=c("a","b","f"))
    peaks2 = GRanges(seqnames = c("1", "2", "2"),
                     IRanges(start = c(1, 100, 200),
                             end = c(1, 101, 201), 
                             names = c("Site1", "Site2", "Site3")),
                     strand="+",feature=c("a","b","f")) #note first
    #chrom diff from above!
    peaks3 = GRanges(seqnames = c("1", "2", "2"),
                     IRanges(start = c(1, 100, 200),
                             end = c(3, 103, 203), 
                             names = c("Site1", "Site2", "Site3")),
                     strand="+",feature=c("a","b","f"))
    peaks4 = GRanges(seqnames = c("1", "2", "2"),
                     IRanges(start = c(1, 100, 200),
                             end = c(2, 101, 201), 
                             names = c("Site1", "Site2", "Site3")),
                     strand="+",feature=c("a","b","f"))
    suppressWarnings(makeVennDiagram(GRangesList(peaks1,peaks2,peaks3),
                                     NameOfPeaks=c("TF1", "TF2","TF3"),
                                     totalTest=100,scaled=F, euler.d=F))
    #is fine, as is below
    suppressWarnings(makeVennDiagram(GRangesList(peaks1,peaks2,peaks3,peaks4),
                                     NameOfPeaks=c("TF1", "TF2","TF3","TF4"),
                                     totalTest=100,scaled=F, euler.d=F))
})