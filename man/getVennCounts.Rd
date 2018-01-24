\name{getVennCounts}
\alias{getVennCounts}
\title{
Obtain Venn Counts for Venn Diagram, internal function for makeVennDigram
}
\description{
Obtain Venn Counts for peak ranges using chromosome ranges or feature field, 
internal function for makeVennDigram
}
\usage{
getVennCounts(..., maxgap = -1L, minoverlap=0L,
    by=c("region", "feature", "base"), 
    ignore.strand=TRUE, connectedPeaks=c("min", "merge", "keepAll"))
}
\arguments{
    \item{\dots}{
    Objects of \link[GenomicRanges:GRanges-class]{GRanges} or 
    \link[IRanges:RangedData-class]{RangedData}: See example below.
}
  \item{maxgap,minoverlap}{
Used in the internal call to \code{findOverlaps()} to detect overlaps.
See \code{?\link[IRanges:findOverlaps-methods]{findOverlaps}} in the \pkg{IRanges} package
for a description of these arguments.
}
\item{by}{region, feature or base, default region. feature means using feature 
field in the RangedData or GRanges for calculating overlap, region means using 
chromosome range for calculating overlap, and base means using calculating 
overlap in nucleotide level.
}
\item{ignore.strand}{
    When set to TRUE, the strand information is ignored in the overlap 
    calculations.
  }
\item{connectedPeaks}{
  If multiple peaks involved in overlapping in several groups, set it to 
  "merge" will count it as only 1, while set it to "min" will count it as 
  the minimal involved peaks in any concered groups 
}
}

\value{
\item{vennCounts}{vennCounts objects containing counts for Venn Diagram 
generation, see details in limma package vennCounts}
}
\author{
Jianhong Ou
}
\seealso{
  \link{makeVennDiagram}, \link{findOverlappingPeaks}
}
\examples{
if(Sys.getenv("USER")=="jianhongou"){
peaks1 = RangedData(IRanges(start = c(967654, 2010897, 2496704), 
                            end = c(967754, 2010997, 2496804), 
                            names = c("Site1", "Site2", "Site3")), 
                   space = c("1", "2", "3"), 
                   strand=as.integer(1), 
                   feature=c("a","b", "c"))
  peaks2 = 
      RangedData(IRanges(start=c(967659, 2010898, 2496700, 3075866, 3123260), 
                         end=c(967869, 2011108, 2496920, 3076166, 3123470), 
                         names = c("t1", "t2", "t3", "t4", "t5")), 
                    space = c("1", "2", "3", "1", "2"), 
                    strand = c(1, 1, -1,-1,1), 
                    feature=c("a","c","d","e", "a"))
    getVennCounts(peaks1,peaks2)
    getVennCounts(peaks1,peaks2, by="feature")
    getVennCounts(peaks1, peaks2, by="base")
}
}
\keyword{ misc }
