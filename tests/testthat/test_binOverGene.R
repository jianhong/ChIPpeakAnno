test_that("binOverGene works not correct", {
  ## by simulation data
  if(Sys.getenv("USER")[1]=="jianhongou"){
  peaks <- toGRanges(TxDb.Hsapiens.UCSC.hg38.knownGene)
  peaks1 <- peaks <- peaks[seqnames(peaks)=="chr1"]
  st <- as.character(strand(peaks))=="+"
  start(peaks)[st] <- start(peaks)[st]-500
  end(peaks)[st] <- start(peaks)[st]
  end(peaks)[!st] <- end(peaks)[!st]+100
  start(peaks)[!st] <- end(peaks)[!st]
  peaks$score <- 1
  start(peaks1)[!st] <- start(peaks1)[!st]-500
  end(peaks1)[!st] <- start(peaks1)[!st]
  end(peaks1)[st] <- end(peaks1)[st]+100
  start(peaks1)[st] <- end(peaks1)[st]
  peaks1$score <- 1
  cvglists <- list(x=coverage(peaks), y=coverage(peaks1))
  bin <- binOverGene(cvglists, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, 
                     upstream.cutoff = 1000)
  expect_true(all(c(11, 19) %in% order(bin$upstream[, "x"], decreasing = TRUE)[1:2]))
  expect_true(all(c(2, 10) %in% order(bin$downstream[, "y"], decreasing = TRUE)[1:2]))
  }
})