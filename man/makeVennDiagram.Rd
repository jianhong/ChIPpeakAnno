\name{makeVennDiagram}
\alias{makeVennDiagram}
\title{
Make Venn Diagram from a list of peaks
}
\description{
Make Venn Diagram from two or more peak ranges, Also calculate p-value to determine whether those peaks overlap significantly.
}
\usage{
makeVennDiagram(Peaks, NameOfPeaks, maxgap = 0L, minoverlap = 1L,
                 totalTest, by = c("region", "feature", "base"),
                 ignore.strand = TRUE, connectedPeaks = c("min",
                 "merge", "keepAll"), method = c("hyperG",
                 "permutation"), TxDb, ...)
}
\arguments{
  \item{Peaks}{
A list of peaks in \link[GenomicRanges:GRanges]{GRanges} format: See example below.
}
  \item{NameOfPeaks}{
Character vector to specify the name of Peaks, e.g., c("TF1", "TF2"). This will be used as label in the Venn Diagram.
}
 \item{maxgap}{
 Non-negative integer. Intervals with a separation of maxgap or less are considered to be overlapping.
}
\item{minoverlap}{
Non-negative integer. Intervals with an overlapping of minoverlap or more are considered to be overlapping.
}
\item{totalTest}{
Numeric value to specify the total number of tests performed to obtain the list of peaks. It should be much larger than the number of peaks in the largest peak set.
}
\item{by}{"region, "feature" or "base", default = "region". feature means using feature field in the GRanges for calculating overlap, region means using chromosome range for calculating overlap, and base means calculating overlap in nucleotide level.
}
\item{ignore.strand}{
    Logical: when set to TRUE, the strand information is ignored in the overlap calculations.
  }
\item{connectedPeaks}{
  If multiple peaks involved in overlapping in several groups, set it to "merge" will count it as only 1, while set it to "min" will count it as the minimal involved peaks in any connected peak group. 
}
\item{method}{method used for p value calculation. hyperG means hypergeometric test and permutation means \link{peakPermTest}}
\item{TxDb}{An object of \link[GenomicFeatures:TxDb]{TxDb}}
\item{\dots}{Additional arguments to be passed to \link[VennDiagram:venn.diagram]{venn.diagram}}
}
\details{
 For customized graph options, please see venn.diagram in VennDiagram package.
}
\value{
In addition to a Venn Diagram produced, a p.value is calculated by hypergeometric test to determine whether the peaks or features are overlapped significantly.
}
\author{
Lihua Julie Zhu, Jianhong Ou
}

\seealso{
 \link{findOverlapsOfPeaks}, \link[VennDiagram:venn.diagram]{venn.diagram}, \link{peakPermTest}
}
\examples{
if (interactive()){
    peaks1 <- GRanges(seqnames=c("1", "2", "3"),
                      IRanges(start=c(967654, 2010897, 2496704),
                              end=c(967754, 2010997, 2496804), 
                              names=c("Site1", "Site2", "Site3")),
                      strand="+",
                      feature=c("a","b","f"))
    peaks2 = RangedData(seqnames=c("1", "2", "3", "1", "2"), 
                        IRanges(start = c(967659, 2010898,2496700,
                                          3075866,3123260),
                                end = c(967869, 2011108, 2496920, 
                                        3076166, 3123470),
                                names = c("t1", "t2", "t3", "t4", "t5")), 
                        strand = c("+", "+", "-", "-", "+"), 
                        feature=c("a","b","c","d","a"))
    makeVennDiagram(list(peaks1, peaks2), NameOfPeaks=c("TF1", "TF2"),
                    totalTest=100,scaled=FALSE, euler.d=FALSE)
    
    makeVennDiagram(list(peaks1, peaks2), NameOfPeaks=c("TF1", "TF2"),
                    totalTest=100)
    
    ###### 4-way diagram using annotated feature instead of chromosome ranges
    
    makeVennDiagram(list(peaks1, peaks2, peaks1, peaks2), 
                    NameOfPeaks=c("TF1", "TF2","TF3", "TF4"), 
                    totalTest=100, by="feature",
                    main = "Venn Diagram for 4 peak lists",
                    fill=c(1,2,3,4))
}
}
\keyword{ graph }
