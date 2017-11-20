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
  }
  \item{upstream, downstream}{Set the feature.gr to upstream and dwonstream 
  from the center of the feature.gr if they are set.}
  \item{n.tile}{The number of tiles to generate for each element of feature.gr, 
    default is 100}
  \item{...}{Not used.}
}
\value{
  A list of matrix. In each matrix, each row record the signals for 
  corresponding feature. rownames of the matrix show the seqnames and coordinates.
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
  featureAlignedSignal(cvglists, feature.gr)
}
\keyword{ misc }
