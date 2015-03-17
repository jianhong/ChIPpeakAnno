test_getAllPeakSequence <- function(){
    #Fixme, Zv9_scaffold.
    	peaks <- RangedData(IRanges(start=c(10,300,60348387),end=c(50,330,60348389), names=c("p1","p2","p0")), space=c("Zv9_scaffold3564","MT","1"))
    	suppressWarnings(seq <- getAllPeakSequence(peaks, upstream = 20, downstream = 20, genome = Drerio))
    	checkEquals(as.character(seq["p0",]$sequence),getSeq(Drerio, GRanges("chr1",IRanges(60348367, 60348388)),as.character=TRUE))
    	checkEquals(as.character(seq["p1",]$sequence),suppressWarnings(getSeq(Drerio, GRanges("Zv9_scaffold3564",IRanges(1, 70)),as.character=TRUE)))
    	checkEquals(as.character(seq["p2",]$sequence),getSeq(Drerio, GRanges("chrM",IRanges(280, 350)),as.character=TRUE))
}