test_that("assignChromosomeRegion works not correct", {
    ##algorithm
    txdb_file <- system.file("extdata", "Biomart_Ensembl_sample.sqlite",
                             package="GenomicFeatures")
    TxDb <- loadDb(txdb_file)
    exons <- exons(TxDb, columns=NULL)
    fiveUTRs <- unique(unlist(fiveUTRsByTranscript(TxDb)))
    Feature.distribution <- 
        assignChromosomeRegion(exons, nucleotideLevel=TRUE, TxDb=TxDb)
    expect_equal(as.integer(Feature.distribution$percentage["Exons"]), 100)
    Feature.distribution <- 
        assignChromosomeRegion(fiveUTRs, nucleotideLevel=FALSE, TxDb=TxDb)
    expect_equal(as.integer(Feature.distribution$percentage["Exons"]), 100)
    expect_equal(as.integer(Feature.distribution$percentage["fiveUTRs"]), 100)
    Feature.distribution <- 
        assignChromosomeRegion(fiveUTRs, precedence="fiveUTRs", TxDb=TxDb)
    expect_equal(as.integer(Feature.distribution$percentage["Exons"]), 0)
    expect_equal(as.integer(Feature.distribution$percentage["fiveUTRs"]), 100)
})