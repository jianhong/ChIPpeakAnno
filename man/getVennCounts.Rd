\name{getVennCounts}
\alias{getVennCounts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain Venn Counts for Venn Diagram, internal function for makeVennDigram
}
\description{
Obtain Venn Counts for peak ranges using chromosome ranges or feature field, internal function for makeVennDigram
}
\usage{
getVennCounts(..., maxgap = 0L, minoverlap=1L, by=c("region", "feature", "base"), 
    ignore.strand=TRUE, connectedPeaks=c("min", "merge", "keepAll"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{\dots}{
    Objects of \link[GenomicRanges:GRanges]{GRanges} or \link[IRanges:RangedData]{RangedData}: See example below.
}
  \item{maxgap}{
Non-negative integer. Intervals with a separation of maxgap or less are considered to be overlapping.
}
\item{minoverlap}{
Non-negative integer. Intervals with an overlapping of minoverlap or more are considered to be overlapping.
}
\item{by}{region, feature or base, default region. feature means using feature field in the RangedData or GRanges for calculating overlap, region means using chromosome range for calculating overlap, and base means using calculating overlap in nucleotide level.
}
\item{ignore.strand}{
    When set to TRUE, the strand information is ignored in the overlap calculations.
  }
\item{connectedPeaks}{
  If multiple peaks involved in overlapping in several groups, set it to "merge" will count it as only 1, while set it to "min" will count it as the minimal involved peaks in any concered groups 
}
}

\value{
\item{vennCounts}{vennCounts objects containing counts for Venn Diagram generation, see details in limma package vennCounts}
}
\author{
Jianhong Ou
}
\seealso{
  \link{makeVennDiagram}, \link{findOverlappingPeaks}
}
\examples{
if(interactive()){
peaks1 = RangedData(IRanges(start = c(967654, 2010897, 2496704), 
                            end = c(967754, 2010997, 2496804), 
                            names = c("Site1", "Site2", "Site3")), 
                   space = c("1", "2", "3"), 
                   strand=as.integer(1), 
                   feature=c("a","b", "c"))
  peaks2 = RangedData(IRanges(start=c(967659, 2010898, 2496700, 3075866, 3123260), 
                              end=c(967869, 2011108, 2496920, 3076166, 3123470), 
                              names = c("t1", "t2", "t3", "t4", "t5")), 
                    space = c("1", "2", "3", "1", "2"), 
                    strand = c(1, 1, -1,-1,1), 
                    feature=c("a","c","d","e", "a"))
    getVennCounts(peaks1,peaks2, maxgap=0)
	getVennCounts(peaks1,peaks2, maxgap=0, by="feature")
    getVennCounts(peaks1, peaks2, maxgap=0, by="base")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
