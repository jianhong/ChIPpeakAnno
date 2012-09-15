\name{makeVennDiagram}
\alias{makeVennDiagram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make Venn Diagram from two peak ranges
}
\description{
Make Venn Diagram from two peak ranges and also calculate p-value for determining whether two peak ranges overlap significantly.
}
\usage{
makeVennDiagram(Peaks, NameOfPeaks, maxgap=0L, minoverlap=1L, totalTest, useFeature=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Peaks}{
RangedDataList: See example below.
}
  \item{NameOfPeaks}{
Character vector to specify the name of Peaks, e.g., c("TF1", "TF2"), this will be used as label in the Venn Diagram.
}
 \item{maxgap}{
 Non-negative integer. Intervals with a separation of maxgap or less are considered to be overlapping.
}
\item{minoverlap}{
Non-negative integer. Intervals with an overlapping of minoverlap or more are considered to be overlapping.}
\item{totalTest}{
Numeric value to specify the total number of tests performed to obtain the list of peaks. It should be much larger than the number of peaks in the largest peak set.
}
\item{useFeature}{TRUE or FALSE, default FALSE, true means using feature field in the RangedData for calculating overlap, false means using chromosome range for calculating overlap.
}
\item{\dots}{Additional arguments to be passed to \link[VennDiagram:venn.diagram]{venn.diagram}}
}
\details{
 For customized graph options, please see venn.diagram in VennDiagram package.
}
\value{
In addition to a Venn Diagram produced, p.value is obtained from hypergeometric test for determining whether the two peak ranges or features overlap significantly.
}
\references{
}
\author{
Lihua Julie Zhu, Jianhong Ou
}
\note{

}

\seealso{
 \link[ChIPpeakAnno:findOverlappingPeaks]{findOverlappingPeaks}, \link[VennDiagram:venn.diagram]{venn.diagram}
}
\examples{
if (interactive())
{
peaks1 = RangedData(IRanges(start = c(967654, 2010897, 2496704),
	 end = c(967754, 2010997, 2496804), names = c("Site1", "Site2", "Site3")),
	 space = c("1", "2", "3"), strand=as.integer(1),feature=c("a","b","f"))
peaks2 = RangedData(IRanges(start = c(967659, 2010898,2496700,3075866,3123260),
	end = c(967869, 2011108, 2496920, 3076166, 3123470),
	names = c("t1", "t2", "t3", "t4", "t5")), 
	space = c("1", "2", "3", "1", "2"), strand = c(1, 1, -1,-1,1), feature=c("a","b","c","d","a"))
makeVennDiagram(RangedDataList(peaks1,peaks2), NameOfPeaks=c("TF1", "TF2"),
 totalTest=100,scaled=F, euler.d=F)

makeVennDiagram(RangedDataList(peaks1,peaks2), NameOfPeaks=c("TF1", "TF2"),
  totalTest=100,useFeature=FALSE)

###### 4-way diagram using annotated feature instead of chromosome ranges

	makeVennDiagram(RangedDataList(peaks1,peaks2, peaks1, peaks2), NameOfPeaks=c("TF1", "TF2","TF3", "TF4"), totalTest=100,useFeature=TRUE,
main = "Venn Diagram for 4 peak lists",fill=c(1,2,3,4))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graph }
