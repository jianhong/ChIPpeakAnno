\name{annoPeaks}
\alias{annoPeaks}
\title{
    Annotate peaks
}
\description{
    Annotate peaks by annoGR object in the given range.
}
\usage{
annoPeaks (peaks, annoData, 
          bindingType=c("bothSidesNSS",
                        "startSite", "endSite", "fullRange", 
                        "bothSidesNearest"), 
          bindingRegion=c(-5000, 500), 
          ignore.peak.strand=TRUE,
          ...)
}
\arguments{
  \item{peaks}{
    peak list, \link[GenomicRanges]{GRanges} object
  }
  \item{annoData}{
    annotation data, \link{annoGR} or \link[GenomicRanges]{GRanges} object
  }
  \item{bindingType}{
    binding relative position. startSite, start position of the feature; endSite,
    end position of the feature; fullRange, whole range of the feature. 
    The strand is considered for startSite and endSite.
    Or if the bindingType is set to bothSidesNSS or bothSidesNearest, 
    annoPeaks will return the nearest start sites or any annotations 
    from both direction of the peak. 
    It will also include all the annotations have overlappings.
  }
  \item{bindingRegion}{
    region range to the annotation. Two integer values. The first one defines 
    the number of nucleotides toward the 5'end and the second value defines the 
    number toward the 3'end, relative to the annotation data. The first one must
    be no bigger than 0. And the second one must be no less than 1. for 
    bindingType equal to startSite or endSite, the downstream site will be 
    restricted in the feature range.
  }
  \item{ignore.peak.strand}{
    ignore the peaks strand or not.
  }
  \item{...}{
    Not used.
  }
}
\value{
    Output is a GRanges object of the annotated peaks.
}
\author{
Jianhong Ou
}
\seealso{
    See Also as \code{\link{annotatePeakInBatch}}, \code{\link{annoGR}}
}
\examples{
    library(EnsDb.Hsapiens.v75)
    data("myPeakList")
    annoGR <- annoGR(EnsDb.Hsapiens.v75)
    seqlevelsStyle(myPeakList) <- seqlevelsStyle(annoGR)
    annoPeaks(myPeakList, annoGR)
}
\keyword{ misc }
