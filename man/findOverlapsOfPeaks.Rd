\name{findOverlapsOfPeaks}
\alias{findOverlapsOfPeaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find the overlapping peaks for two or more peak ranges.
}
\description{
Find the overlapping peaks for two or more (less than five) peak ranges.
}
\usage{
findOverlapsOfPeaks(\dots, maxgap=0L, minoverlap=1L, 
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
\item{ignore.strand}{
    When set to TRUE, the strand information is ignored in the overlap calculations.
  }
\item{connectedPeaks}{
  If multiple peaks involved in overlapping in several groups, set it to "merge" will count it as only 1, while set it to "min" will count it as the minimal involved peaks in any concered groups 
}
}
\details{
Efficiently perform overlap queries with an interval tree implemented in GRanges.
}
\value{
  return value is An object of overlappingPeaks. 
  \item{venn_cnt}{an object of VennCounts}
  \item{peaklist}{a list consists of all overlapping peaks or unique peaks}
}
\references{
1.Interval tree algorithm from: Cormen, Thomas H.; Leiserson, Charles E.; Rivest, Ronald L.; Stein, Clifford. Introduction to Algorithms, second edition, MIT Press and McGraw-Hill. ISBN 0-262-53196-8
2.Zhu L.J. et al. (2010) ChIPpeakAnno: a Bioconductor package to annotate ChIP-seq and ChIP-chip data. BMC Bioinformatics 2010, 11:237doi:10.1186/1471-2105-11-237
}
\author{
Jianhong Ou
}
\seealso{
  \link{annotatePeakInBatch}, \link{makeVennDiagram}, \link{getVennCounts}, \link{findOverlappingPeaks}
}
\examples{
peaks1 <- GRanges(seqnames=c(6,6,6,6,5),
                 IRanges(start=c(1543200,1557200,1563000,1569800,167889600),
                         end=c(1555199,1560599,1565199,1573799,167893599),
                         names=c("p1","p2","p3","p4","p5")),
                 strand="+")
peaks2 <- GRanges(seqnames=c(6,6,6,6,5),
                  IRanges(start=c(1549800,1554400,1565000,1569400,167888600),
                          end=c(1550599,1560799,1565399,1571199,167888999),
                          names=c("f1","f2","f3","f4","f5")),
                  strand="+")
t1 <- findOverlapsOfPeaks(peaks1, peaks2, maxgap=1000)
makeVennDiagram(t1)
t1$venn_cnt
t1$peaklist
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
