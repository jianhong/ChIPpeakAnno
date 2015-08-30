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
findVennCounts(Peaks, NameOfPeaks, maxgap = 0L, minoverlap = 1L, 
               totalTest, useFeature=FALSE)
}
\arguments{
  \item{Peaks}{
RangedDataList: See example below.}
  \item{NameOfPeaks}{
Character vector to specify the name of Peaks, e.g., c("TF1", "TF2"), 
this will be used as label in the Venn Diagram.
}
  \item{maxgap}{
 Non-negative integer. Intervals with a separation of maxgap or less are 
 considered to be overlapping.
}
\item{minoverlap}{
Non-negative integer. Intervals with an overlapping of minoverlap or more are 
considered to be overlapping.}
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
