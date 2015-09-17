\name{featureAlignedDistribution}
\alias{featureAlignedDistribution}
\title{
    plot distribution in given ranges
}
\description{
    plot distribution in the given feature ranges
}
\usage{
featureAlignedDistribution(cvglists, feature.gr, 
                           upstream, downstream, 
                           n.tile=100, zeroAt, ...)
}
\arguments{
  \item{cvglists}{
    Output of \link{featureAlignedSignal} or a 
    list of \link[IRanges]{SimpleRleList} or \link[IRanges]{RleList}
}
  \item{feature.gr}{
    An object of \link[GenomicRanges]{GRanges} with identical width. 
    If the width equal to 1, you can use upstream and downstream to set
    the range for plot.
    If the width not equal to 1, you can use zeroAt to set the zero point
    of the heatmap.
  }
  \item{upstream, downstream}{upstream or dwonstream from the feature.gr.}
  \item{zeroAt}{zero point position of feature.gr}
  \item{n.tile}{The number of tiles to generate for each element of feature.gr, 
    default is 100}
  \item{...}{any paramters could be used by \link[graphics]{matplot}}
}
\value{
  invisible matrix of the plot.
}
\author{
  Jianhong Ou
}
\seealso{
  See Also as \link{featureAlignedSignal}, \link{featureAlignedHeatmap}
}
\examples{
  cvglists <- list(A=RleList(chr1=Rle(sample.int(5000, 100), 
                                      sample.int(300, 100))), 
                   B=RleList(chr1=Rle(sample.int(5000, 100), 
                                      sample.int(300, 100))))
  feature.gr <- GRanges("chr1", IRanges(seq(1, 4900, 100), width=100))
  featureAlignedDistribution(cvglists, feature.gr, zeroAt=50, type="l")
}
\keyword{ misc }
