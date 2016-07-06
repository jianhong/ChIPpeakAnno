\name{featureAlignedExtendSignal}
\alias{featureAlignedExtendSignal}
\title{
    extract signals in given ranges from bam files
}
\description{
    extract signals in the given feature ranges from bam files (DNAseq only). The reads will be 
    extended to estimated fragement length.
}
\usage{
featureAlignedExtendSignal(bamfiles, index = bamfiles, feature.gr, 
                           upstream, downstream, n.tile = 100, 
                           fragmentLength, librarySize, ...)
}
\arguments{
  \item{bamfiles}{
    The file names of the 'BAM' ('SAM' for asBam) files to be processed.
}
  \item{index}{
    The names of the index file of the 'BAM' file being processed; 
    this is given without the '.bai' extension.
}
  \item{feature.gr}{
    An object of \link[GenomicRanges]{GRanges} with identical width.
}
  \item{upstream, downstream}{upstream or dwonstream from the feature.gr.}
  \item{n.tile}{The number of tiles to generate for each element of feature.gr, 
    default is 100}
  \item{fragmentLength}{
    Estimated fragment length.
}
  \item{librarySize}{
    Estimated library size.
}
  \item{\dots}{
   Not used.
}
}
\value{
A list of matrix. In each matrix, each row record the signals for 
  corresponding feature.
}
\author{
Jianhong Ou
}
\seealso{
See Also as \code{\link{featureAlignedSignal}}, \code{\link{estLibSize}}, \code{\link{estFragmentLength}}
}
\examples{
    if(interactive()){
    library(MMDiffBamSubset)
    WT.AB2 <- WT.AB2()
    Null.AB2 <- Null.AB2()
    Resc.AB2 <- Resc.AB2()
    peaks <- WT.AB2.Peaks()
    feature.gr <- toGRanges(peaks, format="MACS")
    feature.gr <- feature.gr[feature.gr$seqnames=="chr1" & 
                             start(feature.gr)>3000000 & 
                             end(feature.gr)<75000000]
    sig <- featureAlignedExtendSignal(c(WT.AB2, Null.AB2, Resc.AB2), 
                               feature.gr=reCenterPeaks(feature.gr, width=1), 
                               upstream = 505,
                               downstream = 505,
                               n.tile=101, 
                               fragmentLength=250,
                               librarySize=1e9)
    featureAlignedHeatmap(sig, reCenterPeaks(feature.gr, width=1010), 
                          zeroAt=.5, n.tile=101)
    }
}
\keyword{ misc }
