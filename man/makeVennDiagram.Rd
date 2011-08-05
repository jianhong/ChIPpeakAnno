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
makeVennDiagram(Peaks, NameOfPeaks, maxgap=0, totalTest, cex = 1.5,
 counts.col = "red", useFeature=FALSE)
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
\item{totalTest}{
Numeric value to specify the total number of tests performed to obtain the list of peaks. It should be much larger than the number of peaks in the largest peak set.
}
  \item{cex}{
Numerical value giving the amount by which the contrast names should be scaled on the plot relative to the default.plotting text. See par.}
  \item{counts.col}{
 optional vector of color specifications defining the colors by which the circles should be drawn. See par.
}
\item{useFeature}{TRUE or FALSE, default FALSE, true means using feature field in the RangedData for calculating overlap, false means using chromosome range for calculating overlap.
}
}
\details{

}
\value{
In addition to a Venn Diagram produced, p.value is obtained from hypergeometric test for determining whether the two peak ranges or fetures overlap significantly.
}
\references{
}
\author{
Lihua Julie Zhu
}
\note{

}

\seealso{
 findOverlappingPeaks
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
 totalTest=100)

makeVennDiagram(RangedDataList(peaks1,peaks2), NameOfPeaks=c("TF1", "TF2"),
  totalTest=100,useFeature=FALSE)

###### 4-way diagram using annotated feature instead of chromosome ranges

	makeVennDiagram(RangedDataList(peaks1,peaks2, peaks1, peaks2), NameOfPeaks=c("TF1", "TF2","TF1", "TF2"),
  totalTest=100,useFeature=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graph }
