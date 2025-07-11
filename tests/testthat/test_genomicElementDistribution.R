test_that("assignChromosomeRegion works not correct", {
    ##algorithm
    txdb_file <- system.file("extdata", "Biomart_Ensembl_sample.sqlite",
                             package="GenomicFeatures")
    TxDb <- loadDb(txdb_file)
    exons <- exons(TxDb, columns=NULL)
    fiveUTRs <- unique(unlist(fiveUTRsByTranscript(TxDb)))
    Feature.distribution <- 
        genomicElementDistribution(exons, nucleotideLevel=TRUE, TxDb=TxDb)
    expect_true(all(Feature.distribution$peaks$ExonIntron=='exon'))
    Feature.distribution <- 
      genomicElementDistribution(exons, TxDb=TxDb)
    Feature.distribution <- 
      genomicElementDistribution(fiveUTRs, TxDb=TxDb)
    expect_true(all(Feature.distribution$peaks$Exons=='utr5'))
    seqlevelsStyle(fiveUTRs) <- 'UCSC'
    Feature.distribution <- 
      genomicElementDistribution(fiveUTRs, TxDb=TxDb)
    expect_equal(seqlevelsStyle(Feature.distribution$peaks), 'UCSC')
    expect_true(all(Feature.distribution$peaks$Exons=='utr5'))
})