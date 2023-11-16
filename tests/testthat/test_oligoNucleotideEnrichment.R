test_that("oligoNucleotideEnrichment works not correct", {
filepath =system.file("extdata", "examplePattern.fa", package="ChIPpeakAnno")
peaks = GRanges(seqnames=c("chr17","chr3"),
                 IRanges(start=c(41275784,10076141),
                         end=c(41276382,10076732),
                         names=paste0("peak",1:2)))

result <- oligoNucleotideEnrichment(filepath=filepath,
peaks=peaks,
genome=Hsapiens,
methodBackground="selectChromRandomly",
times=20,
alpha=0.1)

expect_equal(result[1:2,4], c(0.6869590,0.7573282),tolerance=1e-5)
expect_equal(result[1:2,3], c(0.0009829383,0.0033579623),tolerance=1e-6)
expect_equal(result[1:2,2], c(1181,1181))
expect_equal(result[1:2,1], c(1,3))

})