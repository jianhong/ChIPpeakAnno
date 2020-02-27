test_that("binOverRegions or binOverGene works not correct", {
    ##algorithm
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    if(Sys.getenv("USER")[1]=="jianhongou"){
        TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
        x <- isActiveSeq(TxDb)
        x[seq_along(x)] <- FALSE
        x["chrX"] <- TRUE
        isActiveSeq(TxDb) <- x
        utr5 <- fiveUTRsByTranscript(TxDb)
        utr3 <- threeUTRsByTranscript(TxDb)
        CDS <- cds(TxDb)
        gene <- genes(TxDb)
        exon <- exons(TxDb)
        txs <- transcripts(TxDb)
        utr5 <- unlist(utr5)
        utr3 <- unlist(utr3)
        prom <- promoters(gene, upstream=1000, downstream=0)
        anno <- list("5UTR"=utr5,
                     "3UTR"=utr3,
                     "CDS" =CDS,
                     "gene"=gene,
                     "exon"=exon,
                     "txs" =txs,
                     "prom"=prom)
        anno.rd <- lapply(anno, reduce)
        cvglists.rd <- lapply(anno.rd, coverage)
        d1 <- binOverRegions(cvglists.rd, TxDb)
        d2 <- binOverGene(cvglists.rd, TxDb)
        plotBinOverRegions(d1, main="binOverRegions.reduced")
        plotBinOverRegions(d2, main="binOverGene.reduced")
    }
})