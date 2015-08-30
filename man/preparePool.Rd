\name{preparePool}
\alias{preparePool}
\title{
    prepare data for permutation test
}
\description{
    prepare data for permutation test \link{peakPermTest}
}
\usage{
preparePool(TxDb, template, bindingDistribution, 
            bindingType = c("TSS", "geneEnd"), 
            featureType = c("transcript", "exon"), 
            seqn = NA)
}
\arguments{
  \item{TxDb}{
an object of \link[GenomicFeatures]{TxDb}
}
  \item{template}{
an object of \link[GenomicRanges]{GRanges}
}
  \item{bindingDistribution}{
an object of \link{bindist}
}
  \item{bindingType}{
the relevant position to features
}
  \item{featureType}{
feature type, transcript or exon.
}
  \item{seqn}{
seqnames. If given, the pool for permutation will be restrict in the 
given chromosomes.
}
}
\value{
a list with two elements,
 grs, a list of \link[GenomicRanges]{GRanges}.
 N, the numbers of elements should be drawn from in each GRanges.
}
\seealso{
    \link{peakPermTest}, \link{bindist}
}
\author{
Jianhong Ou
}
\examples{
    if(interactive()){
        path <- system.file("extdata", package="ChIPpeakAnno")
        peaksA <- toGRanges(file.path(path, "peaks.narrowPeak"), 
                            format="narrowPeak")
        peaksB <- toGRanges(file.path(path, "MACS2_peaks.xls"), format="MACS2")
        library(TxDb.Hsapiens.UCSC.hg19.knownGene)
        ppp <- preparePool(TxDb.Hsapiens.UCSC.hg19.knownGene, 
                           peaksA, bindingType="TSS",
                           featureType="transcript")
    }
}

\keyword{ misc }