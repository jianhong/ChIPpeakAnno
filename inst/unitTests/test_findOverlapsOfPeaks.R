test_findOverlapsOfPeaks<-function(){
    ##2 way
    p1 <- GRanges("1", 
                  IRanges(start=c(1, 3, 5, 8, 10), 
                          end=c(2, 4, 6, 9, 11), 
                          names=letters[1:5]), 
                  strand="*")
    p2 <- GRanges("1", 
                  IRanges(start=c(2, 4, 9), 
                          end=c(3, 5, 10), 
                          names=letters[1:3]), 
                  strand="*")
    ol <- findOverlapsOfPeaks(p1, p2, maxgap=0)
    checkEquals(ol$venn_cnt[,"Counts"], c(0,0,0,3))
    ##3 way
    p3 <- GRanges("1", 
                  IRanges(start=c(3), end=c(4), names=letters[1]), 
                  strand="*")
    ol <- findOverlapsOfPeaks(p1, p2, p3, maxgap=0)
    checkEquals(ol$venn_cnt[,"Counts"], c(0,0,0,0,0,0,1,1))
    ##4 way
    p4 <- GRanges("1", 
                  IRanges(start=c(13), end=c(14), names=letters[1]), 
                  strand="*")
    ol <- findOverlapsOfPeaks(p1, p2, p3, p4, maxgap=0)
    checkEquals(ol$venn_cnt[,"Counts"], 
                c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0))
    ##5 way
    p5 <- GRanges("1", 
                  IRanges(start=c(14), end=c(15), names=letters[1]), 
                  strand="*")
    ol <- findOverlapsOfPeaks(p1, p2, p3, p4, p5, maxgap=0)
    ## minoverlap > 2
    p1 <- GRanges(rep(1:8, each=3), 
                  IRanges(start=rep(c(1, 5, 8), 8), end=rep(c(4, 7, 11), 8), 
                          names=letters[1:24]), 
                  strand="*", id=LETTERS[1:24])
    p2 <- GRanges(rep(1:8, each=2), 
                  IRanges(start=rep(c(2, 6), 8), end=rep(c(5, 10), 8), 
                          names=letters[1:16]), 
                  strand="*", id=LETTERS[1:16])
    ol <- findOverlapsOfPeaks(p1, p2, maxgap=0, minoverlap=3)
    checkEquals(p1[gsub("p1__", "", 
                        unlist(ol$peaklist$p1$peakNames))]$id, 
                LETTERS[1:8*3-1])
    d <- ol$peaklist$`p1///p2`$peakNames
    d <- do.call(rbind, lapply(d, function(.ele) .ele[order(.ele)]))
    d <- cbind(p1[gsub("p1__", "", d[,1])]$id, p2[gsub("p2__", "", d[,2])]$id)
    d <- apply(d, 1, sort)
    checkEquals(paste(d[1,], d[2,]), 
                paste(LETTERS[1:16], LETTERS[1:24][-(1:8*3-1)]))
}