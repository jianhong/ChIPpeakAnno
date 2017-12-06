\name{bdp}
\alias{bdp}
\title{
    obtain the peaks near bi-directional promoters
}
\description{
    Obtain the peaks near bi-directional promoters. Also output percent of 
  peaks near bi-directional promoters.
}
\usage{
bdp (peaks, annoData, maxgap=2000L, ...)
}
\arguments{
  \item{peaks}{
    peak list, \link[GenomicRanges]{GRanges} object
  }
  \item{annoData}{
    annotation data, \link{annoGR} object
  }
  \item{maxgap}{
    maxgap between peak and TSS
  }
  \item{...}{
    Not used.
  }
}
\value{
    Output is a list of GRanges object of the peaks near bi-directional promoters.
}
\author{
Jianhong Ou
}
\seealso{
    See Also as \code{\link{annoPeaks}}, \code{\link{annoGR}}
}
\examples{
  if(Sys.getenv("USER")=="jianhongou"){
    library(EnsDb.Hsapiens.v75)
    data("myPeakList")
    annoGR <- annoGR(EnsDb.Hsapiens.v75)
    seqlevelsStyle(myPeakList) <- seqlevelsStyle(annoGR)
    ChIPpeakAnno:::bdp(myPeakList, annoGR)
  }
}
\keyword{ misc }
