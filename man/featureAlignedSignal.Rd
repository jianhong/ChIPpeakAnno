\name{featureAlignedSignal}
\alias{featureAlignedSignal}
\title{
    extract signals in given ranges
}
\description{
    extract signals in the given feature ranges
}
\usage{
featureAlignedSignal(cvglists, feature.gr, 
                     upstream, downstream, 
                     n.tile=100, ...)
}
\arguments{
  \item{cvglists}{
    List of \link[IRanges]{SimpleRleList} or \link[IRanges]{RleList}
  }
  \item{feature.gr}{
    An object of \link[GenomicRanges]{GRanges} with identical width. 
    If the width equal to 1, you can use upstream and downstream to set
    the range for plot.
    If the width not equal to 1, you can use zeroAt to set the zero point
    of the heatmap.
  }
  \item{upstream, downstream}{upstream or dwonstream from the feature.gr.}
  \item{n.tile}{The number of tiles to generate for each element of feature.gr, 
    default is 100}
  \item{...}{Not used.}
}
\value{
  A list of matrix. In each matrix, each row record the signals for 
  corresponding feature.
}
\author{
  Jianhong Ou
}
\seealso{
  See Also as \link{featureAlignedHeatmap}, \link{featureAlignedDistribution}
}
\examples{
  cvglists <- list(A=RleList(chr1=Rle(sample.int(5000, 100), 
                                      sample.int(300, 100))), 
                   B=RleList(chr1=Rle(sample.int(5000, 100), 
                                      sample.int(300, 100))))
  feature.gr <- GRanges("chr1", IRanges(seq(1, 4900, 100), width=100))
  featureAlignedSignal(cvglists, feature.gr, zeroAt=50, type="l")
}
\keyword{ misc }
