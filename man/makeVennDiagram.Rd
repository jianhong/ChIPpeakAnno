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
makeVennDiagram(Peaks1, Peaks2, overlappingPeaks, totalTest, NameOfPeaks1 = "TF1", NameOfPeaks2 = "TF2", cex = 1, counts.col = "red")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Peaks1}{
RangedData: See example below.
}
  \item{Peaks2}{
RangedData: See example below.
}
  \item{overlappingPeaks}{
Data frame: obtained from findOverlappingPeaks(Peaks1, Peaks2, maxgap=100,multiple= FALSE), see example below.
}
  \item{totalTest}{
Numeric value to specify the total number of tests performed to obtain the list of peaks.
}
  \item{NameOfPeaks1}{
Name of Peaks1, e.g., TF1, this will be used as label in the Venn Diagram.
}
  \item{NameOfPeaks2}{
Name of Peaks2, e.g., TF2, this will be used as label in the Venn Diagram.
}
  \item{cex}{
Numerical value giving the amount by which the contrast names should be scaled on the plot relative to the default.plotting text. See par.}
  \item{counts.col}{
 optional vector of color specifications defining the colors by which the circles should be drawn. See par.
}
}
\details{
This is a wrapper function for vennDiagram from limma package.
}
\value{
In addition to a Venn Diagram produced, p.value is obtained from hypergeometric test for determining whether the two peak ranges overlap significantly.
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
	myexp = RangedData(IRanges(start = c(967654, 2010897, 2496704), end = c(967754, 2010997, 2496804), names = c("Site1", "Site2", "Site3")), space = c("1", "2", "3"))	literature = RangedData(IRanges(start = c(967659, 2010898,  2496700, 3075866, 3123260), end = c(967869, 2011108, 2496920, 3076166, 3123470), names = c("t1", "t2", "t3", "t4", "t5")), space = c("1", "2", "3", "1", "2"), strand = c(1, 1, -1,-1,1))
	t1 =findOverlappingPeaks(myexp, literature, maxgap=100, multiple=F) 
	r = t1$OverlappingPeaks
	makeVennDiagram(myexp, literature, r, totalTest=100)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graph }
