test_that("getVennCounts works not correct", {
    ##2 way
    p1 <- GRanges("1", 
                  IRanges(start=c(1, 3, 5, 8, 10), 
                          end=c(2, 4, 6, 9, 11), 
                          names=letters[1:5]), 
                  strand="*", 
                  feature=c("a","c","d","e", "a"))
    p2 <- GRanges("1", 
                  IRanges(start=c(2, 4, 9), 
                          end=c(3, 5, 10), 
                          names=letters[1:3]), 
                  strand="*", 
                  feature=c("a","b", "c"))
    venn_cnt <- getVennCounts(p1, p2, maxgap=0, by="region")
    expect_equal(venn_cnt[,"Counts"], c(0,0,0,3))
    venn_cnt <- getVennCounts(p1, p2, maxgap=0, by="feature")
    venn_cnt <- getVennCounts(p1, p2, maxgap=0, by="base")
    ##3 way
    p3 <- GRanges("1", 
                  IRanges(start=c(3), end=c(4), names=letters[1]), 
                  strand="*", 
                  feature=c("a"))
    venn_cnt <- getVennCounts(p1, p2, p3, maxgap=0, by="region")
    expect_equal(venn_cnt[,"Counts"], c(0,0,0,0,0,0,1,1))
    venn_cnt <- getVennCounts(p1, p2, p3, maxgap=0, by="feature")
    ##4 way
    p4 <- GRanges("1", 
                  IRanges(start=c(13), end=c(14), names=letters[1]), 
                  strand="*", 
                  feature=c("a"))
    venn_cnt <- getVennCounts(p1, p2, p3, p4, maxgap=0, by="region")
    expect_equal(venn_cnt[,"Counts"], 
                c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0))
    venn_cnt <- getVennCounts(p1, p2, p3, p4, maxgap=0, by="feature")
})