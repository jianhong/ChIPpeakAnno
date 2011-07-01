\name{BED2RangedData}
\alias{BED2RangedData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ convert BED format to RangedData }
\description{
  convert BED format to RangedData}
\usage{
BED2RangedData(data.BED,header=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.BED}{ BED format data frame, please refer to http://genome.ucsc.edu/FAQ/FAQformat#format1 for details }
  \item{header}{TRUE or FALSE, default to FALSE, indicates whether data.BED file has BED header} 
}
\details{
  
}
\value{
   RangedData with slot start holding the start position of the feature, slot end holding the end position of the feature, slot names holding the id of the feature, slot space holding the chromosome location where the feature is located. In addition, the following variables are included.
    \item{\code{strand}}{1 for positive strand and -1 for negative strand where the feature is located. Default to 1 if not present in the BED formated data frame}
}
\references{  }
\author{ Lihua Julie Zhu }
\note{ 
For converting the peakList in BED format to RangedData before calling annotatePeakInBatch function
}
\seealso{  }
\examples{
	test.bed = data.frame(cbind(chrom = c("1", "2"), chromStart=c("100", "1000"),
	chromEnd=c("200", "1100"), name=c("peak1", "peak2")))
	test.rangedData = BED2RangedData(test.bed)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }% __ONLY ONE__ keyword per line
