\name{findVennCounts}
\alias{findVennCounts}
\title{
Obtain Venn Counts for Venn Diagram, internal function for makeVennDigram
}
\description{
Obtain Venn Counts for two peak ranges using chromosome ranges or feature 
field, internal function for makeVennDigram
}
\usage{
findVennCounts(Peaks, NameOfPeaks, maxgap = -1L, minoverlap = 0L,
               totalTest, useFeature=FALSE)
}
\arguments{
  \item{Peaks}{
RangedDataList: See example below.}
  \item{NameOfPeaks}{
Character vector to specify the name of Peaks, e.g., c("TF1", "TF2"), 
this will be used as label in the Venn Diagram.
}
  \item{maxgap,minoverlap}{
Used in the internal call to \code{findOverlaps()} to detect overlaps.
See \code{?\link[IRanges:findOverlaps-methods]{findOverlaps}} in the \pkg{IRanges} package
for a description of these arguments.
}
  \item{totalTest}{
Numeric value to specify the total number of tests performed to obtain the 
list of peaks.
}
\item{useFeature}{TRUE or FALSE, default FALSE, true means using feature field 
in the RangedData for calculating overlap, false means using chromosome range 
for calculating overlap.
}
}
\value{
 \item{p.value }{hypergeometric testing result}
\item{vennCounts}{vennCounts objects containing counts for Venn Diagram 
generation, see details in limma package vennCounts}

}
\author{
Lihua Julie Zhu
}
\seealso{
makeVennDiagram
}
\keyword{ misc }
