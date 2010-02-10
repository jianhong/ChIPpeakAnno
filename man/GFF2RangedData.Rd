\name{GFF2RangedData}
\alias{GFF2RangedData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ convert GFF format to RangedData }
\description{
  convert GFF format to RangedData}
\usage{
GFF2RangedData(data.GFF,header=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.GFF}{ GFF format data frame, please refer to http://genome.ucsc.edu/FAQ/FAQformat#format3 for details }
  \item{header}{TRUE or FALSE, default to FALSE, indicates whether data.GFF file has GFF header} 
}
\details{
  
}
\value{
   RangedData with slot start holding the start position of the feature, slot end holding the end position of the feature, slot names holding the id of the feature, slot space holding the chromosome location where the feature is located. In addition, the following variables are included.
    \item{\code{strand}}{1 for positive strand and -1 for negative strand where the feature is located.}
}
\references{  }
\author{ Lihua Julie Zhu }
\note{ 
For converting the peakList in GFF format to RangedData before calling annotatePeakInBatch function
}
\seealso{  }
\examples{
	test.GFF = data.frame(cbind(seqname  = c("chr1", "chr2"), source=rep("Macs", 2), feature=rep("peak", 2), start=c("100", "1000"), end=c("200", "1100"), score=c(60, 26), strand=c(1, -1), frame=c(".", 2), group=c("peak1", "peak2")))
	test.rangedData = GFF2RangedData(test.GFF)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }% __ONLY ONE__ keyword per line
