test_assignChromosomeRegion<-function(){
    ##algorithm
    txdb_file <- system.file("extdata", "Biomart_Ensembl_sample.sqlite",
                             package="GenomicFeatures")
    TxDb <- loadDb(txdb_file)
    exons <- exons(TxDb, columns=NULL)
    fiveUTRs <- unique(unlist(fiveUTRsByTranscript(TxDb)))
    Feature.distribution <- assignChromosomeRegion(exons, nucleotideLevel=TRUE, TxDb=TxDb)
    checkEqualsNumeric(Feature.distribution$percentage["Exons"], 100, tolerance=1e-05)
    Feature.distribution <- assignChromosomeRegion(fiveUTRs, nucleotideLevel=FALSE, TxDb=TxDb)
    checkEqualsNumeric(Feature.distribution$percentage["Exons"], 100, tolerance=1e-05)
    checkEqualsNumeric(Feature.distribution$percentage["fiveUTRs"], 100, tolerance=1e-05)
    Feature.distribution <- assignChromosomeRegion(fiveUTRs, precedence="fiveUTRs", TxDb=TxDb)
    checkEqualsNumeric(Feature.distribution$percentage["Exons"], 0, tolerance=1e-05)
    checkEqualsNumeric(Feature.distribution$percentage["fiveUTRs"], 100, tolerance=1e-05)
}