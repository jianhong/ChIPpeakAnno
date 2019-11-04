\name{ChIPpeakAnno-deprecated}
\alias{ChIPpeakAnno-deprecated}
\title{
Deprecated Functions in Package ChIPpeakAnno
}
\description{
These functions are provided for compatibility with older versions of R only, 
and may be defunct as soon as the next release.
}
\usage{
  findOverlappingPeaks(Peaks1, Peaks2, maxgap = -1L,
                    minoverlap=0L, multiple = c(TRUE, FALSE),
                    NameOfPeaks1 = "TF1", NameOfPeaks2 = "TF2", 
                    select=c("all", "first","last","arbitrary"), 
                    annotate = 0, ignore.strand=TRUE, 
                    connectedPeaks=c("min", "merge"), \dots)
}
\arguments{
  \item{Peaks1}{GRanges: See example below.}
  \item{Peaks2}{GRanges: See example below.}
  \item{maxgap,minoverlap}{
Used in the internal call to \code{findOverlaps()} to detect overlaps.
See \code{?\link[IRanges:findOverlaps-methods]{findOverlaps}} in the \pkg{IRanges} package
for a description of these arguments.
}
  \item{multiple}{TRUE or FALSE: TRUE may return multiple overlapping peaks in 
  Peaks2 for one peak in Peaks1; FALSE will return at most one overlapping 
  peaks in Peaks2 for one peak in Peaks1. This parameter is kept for backward 
  compatibility, please use select.}
  \item{NameOfPeaks1}{Name of the Peaks1, used for generating column name.}
  \item{NameOfPeaks2}{Name of the Peaks2, used for generating column name.}
  \item{select}{all may return multiple overlapping peaks, first will return 
  the first overlapping peak, last will return the last overlapping peak and 
  arbitrary will return one of the overlapping peaks.
}
  \item{annotate}{
Include overlapFeature and shortestDistance in the OverlappingPeaks or not. 
1 means yes and 0 means no. Default to 0.}
\item{ignore.strand}{
    When set to TRUE, the strand information is ignored in the overlap 
    calculations.
  }
\item{connectedPeaks}{
  If multiple peaks involved in overlapping in several groups, set it to 
  "merge" will count it as only 1, while set it to "min" will count it as 
  the minimal involved peaks in any concered groups 
}
\item{\dots}{
    Objects of \link[GenomicRanges:GRanges-class]{GRanges}: See also 
    \code{\link{findOverlapsOfPeaks}}.
}
}
\details{
  findOverlappingPeaks is now deprecated wrappers for 
  \code{\link{findOverlapsOfPeaks}}
}
\seealso{
  \code{\link{Deprecated}}, \link{findOverlapsOfPeaks}, \link{toGRanges}
}
