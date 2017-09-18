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
  BED2RangedData(data.BED,header=FALSE, \dots)
  GFF2RangedData(data.GFF,header=FALSE, \dots)
}
\arguments{
  \item{Peaks1}{RangedData: See example below.}
  \item{Peaks2}{RangedData: See example below.}
  \item{maxgap,minoverlap}{
Used in the internal call to \code{findOverlaps()} to detect overlaps.
See \code{?\link[IRanges]{findOverlaps}} in the \pkg{IRanges} package
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
    Objects of \link[GenomicRanges:GRanges]{GRanges} or 
    \link[IRanges:RangedData]{RangedData}: See also 
    \code{\link{findOverlapsOfPeaks}}. 
    Or any parameter need to be passed into read.delim function for 
    2RangedData function.
}
\item{header}{TRUE or FALSE, default to FALSE, 
indicates whether data file has header}
\item{data.BED}{ BED format data frame or BED filename, 
please refer to http://genome.ucsc.edu/FAQ/FAQformat#format1 for details }
\item{data.GFF}{ GFF format data frame or GFF file name, 
please refer to http://genome.ucsc.edu/FAQ/FAQformat#format3 for details }
}
\details{
  findOverlappingPeaks is now deprecated wrappers for 
  \code{\link{findOverlapsOfPeaks}}
}
\seealso{
  \code{\link{Deprecated}}, \link{findOverlapsOfPeaks}, \link{toGRanges}
}
