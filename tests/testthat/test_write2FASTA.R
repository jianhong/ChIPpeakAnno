test_that("write2FASTA works not correct", {
    ##GRanges
    peaksWithSequences = GRanges("chr1", 
                                 IRanges(start=c(1000, 2000), 
                                         end=c(1010, 2010), 
                                         names=c("id1", "id2")), 
                                 sequence= c("CCCCCCCCGGGGG", 
                                             "TTTTTTTAAAAAA"))
    write2FASTA(peaksWithSequences, tempfile(fileext = "testWrite2FASTA.fa"))
})