\name{findVennCounts}
\alias{findVennCounts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain Venn Counts for Venn Diagram, internal function for makeVennDigram
}
\description{
Obtain Venn Counts for two peak ranges using chromosome ranges or feature field, internal function for makeVennDigram
}
\usage{
findVennCounts(Peaks, NameOfPeaks, maxgap = 0, totalTest, useFeature=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Peaks}{
RangedDataList: See example below.}
  \item{NameOfPeaks}{
Character vector to specify the name of Peaks, e.g., c("TF1", "TF2"), this will be used as label in the Venn Diagram.
}
  \item{maxgap}{
 Non-negative integer. Intervals with a separation of maxgap or less are considered to be overlapping.
}
  \item{totalTest}{
Numeric value to specify the total number of tests performed to obtain the list of peaks.
}
\item{useFeature}{TRUE or FALSE, default FALSE, true means using feature field in the RangedData for calculating overlap, false means using chromosome range for calculating overlap.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 \item{p.value }{hypergeometric testing result}
\item{vennCounts}{vennCounts objects containing counts for Venn Diagram generation, see details in limma package vennCounts}

}
\references{

}
\author{
Lihua Julie Zhu
}
\note{
	if (interactive())
	{
	peaks1 = RangedData(IRanges(start = c(967654, 2010897, 2496704), end = c(967754, 2010997, 2496804), names = c("Site1", "Site2", "Site3")), space = c("1", "2", "3"), strand=as.integer(1), feature=c("a","b", "c"))
	peaks2 = RangedData(IRanges(start = c(967659, 2010898,  2496700, 3075866, 3123260), end = c(967869, 2011108, 2496920, 3076166, 3123470), names = c("t1", "t2", "t3", "t4", "t5")), space = c("1", "2", "3", "1", "2"), strand = c(1, 1, -1,-1,1), feature=c("a","c","d","e", "a"))
	findVennCounts(RangedDataList(peaks1,peaks2), NameOfPeaks=c("TF1", "TF2"), maxgap=0,totalTest=100, useFeature=FALSE)
	findVennCounts(RangedDataList(peaks1,peaks2), NameOfPeaks=c("TF1", "TF2"), maxgap=0,totalTest=100, useFeature=TRUE)
	}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
makeVennDiagram
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
