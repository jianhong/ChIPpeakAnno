\name{reCenterPeaks}
\alias{reCenterPeaks}
\title{
    re-center the peaks
}
\description{
    Create a new list of peaks based on the peak centers of given list.
}
\usage{
reCenterPeaks(peaks, width=2000L, ...)
}
\arguments{
  \item{peaks}{
    An object of \link[GenomicRanges]{GRanges} or \link{annoGR}.
}
  \item{width}{
    The width of new peaks
  }
  \item{...}{Not used.}
}
\value{
    An object of GRanges.
}
\author{
    Jianhong Ou
}

\examples{
    reCenterPeaks(GRanges("chr1", IRanges(1, 10)), width=2)
}

\keyword{ misc }