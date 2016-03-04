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
          bindingType=c("nearestBiDirectionalPromoters",
                        "startSite", "endSite", "fullRange"), 
          bindingRegion=c(-5000, 5000), 
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
    Specifying the criteria to associate peaks with annotation. Here is how to use it together with the parameter bindingRegion
    \itemize{
        \item To obtain peaks within 5kb upstream and up to 3kb downstream of TSS within the gene body, set bindingType = "startSite" and bindingRegion = c(-5000, 3000)
        \item To obtain  peaks  up to 5kb upstream within the gene body and 3kb downstream of gene/Exon End, set bindingType = "endSite" and bindingRegion = c(-5000, 3000)
        \item To obtain peaks from 5kb upstream to 3kb downstream of genes/Exons , set  bindingType = "fullRange" and bindingRegion = c(-5000, 3000)
        \item To obtain peaks with nearest bi-directional promoters within 5kb upstream and 3kb downstream of TSS, set  bindingType = "nearestBiDirectionalPromoters" and bindingRegion = c(-5000, 3000)
    }
    \describe{
        \item{startSite}{start position of the feature (strand is considered)} 
        \item{endSite}{end position of the feature (strand is considered)} 
        \item{fullRange}{whole range of the feature} 
        \item{nearestBiDirectionalPromoters}{nearest promoters from both direction of the peaks (strand is considered). It will report bidirectional promoters if there are promoters in both directions in the given region (defined by bindingRegion). Otherwise, it will report the closest promoter in one direction.}
    }
  }
  \item{bindingRegion}{
    Annotation range used together with bindingType, which is a vector with two integer values, default to c (-5000, 5000). The first one must be no bigger than 0. And the sec ond one must be no less than 1.  For details, see bindingType.
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
