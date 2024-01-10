test_that("summarizePatternInPeaks does not work correctly", {
filepath <- system.file("extdata", "examplePattern.fa", package = "ChIPpeakAnno")

# Check the motif occurrences in one example peak
# Skip checking for p-value calculation
peaks <- GRanges(seqnames = c("chr17"),
                 IRanges(start = c(41275784),
                         end = c(41296382),
                         names = paste0("peak", 1)))
peaks_seq <- getAllPeakSequence(peaks, genome = BSgenome.Hsapiens.UCSC.hg19, 
                                upstream = 0, downstream = 0)$sequence

# Test for each motif separately
patternVec <- readDNAStringSet(filepath, "fasta", use.names = TRUE)

for (i in 1:length(patternVec)) {
  pattern <- patternVec[i]
  temfile <- tempfile()
  writeXStringSet(pattern, temfile)
  
  # Without reverse complement:
  result <- summarizePatternInPeaks(patternFilePath = temfile,
                                    peaks = peaks, 
                                    BSgenomeName = BSgenome.Hsapiens.UCSC.hg19,
                                    revcomp = FALSE)
  
  if (pattern == "AACCNM") {
    res <- gregexec("AACC[A|T|C|G][A|C]", peaks_seq, perl = TRUE)
    expect_equal(length(res[[1]]), nrow(result$motif_occurrence))
  } else if (pattern == "GGNCCK") {
    res <- gregexec("GG[A|T|C|G]CC[G|T]", peaks_seq, perl = TRUE)
    expect_equal(length(res[[1]]), nrow(result$motif_occurrence))
  }
  
  # With reverse complement:
  result <- summarizePatternInPeaks(patternFilePath = temfile,
                                    peaks = peaks, 
                                    BSgenomeName = BSgenome.Hsapiens.UCSC.hg19,
                                    revcomp = TRUE)
  
  if (pattern == "AACCNM") {
    res1 <- gregexec("AACC[A|T|C|G][A|C]", peaks_seq, perl = TRUE)
    res2 <- gregexec("[T|G][G|C|A|T]GGTT", peaks_seq, perl = TRUE)
    expect_equal(length(res1[[1]]) + length(res2[[1]]), nrow(result$motif_occurrence))
  } else if (pattern == "GGNCCK") {
    res1 <- gregexec("GG[A|T|C|G]CC[G|T]", peaks_seq, perl = TRUE)
    res2 <- gregexec("[C|A]GG[C|G|A|T]CC", peaks_seq, perl = TRUE)
    expect_equal(length(res1[[1]]) + length(res2[[1]]), nrow(result$motif_occurrence))
  }
}

})