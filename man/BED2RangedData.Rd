\name{BED2RangedData}
\alias{BED2RangedData}
\alias{BED2RangedData-deprecated}
\title{ Convert BED format to RangedData }
\description{
  Convert BED format to RangedData. This function will be depreciated. }
\usage{
BED2RangedData(data.BED,header=FALSE, ...)
}
\arguments{
  \item{data.BED}{ BED format data frame or BED filename, please refer to 
  http://genome.ucsc.edu/FAQ/FAQformat#format1 for details }
  \item{header}{TRUE or FALSE, default to FALSE, indicates whether data.BED 
  file has BED header}
  \item{...}{any parameter need to be passed into read.delim function}
}
\value{
   RangedData with slot start holding the start position of the feature, slot 
   end holding the end position of the feature, slot names holding the id of 
   the feature, slot space holding the chromosome location where the feature 
   is located. In addition, the following variables are included.
    \item{\code{strand}}{1 for positive strand and -1 for negative strand where 
    the feature is located. Default to 1 if not present in the BED formated 
    data frame}
}
\author{ Lihua Julie Zhu }
\note{ 
For converting the peakList in BED format to RangedData before calling 
annotatePeakInBatch function
}
\seealso{
    See also as \code{\link{toGRanges}}.
}
\examples{
    test.bed = data.frame(cbind(chrom = c("1", "2"), 
                          chromStart=c("100", "1000"),
                          chromEnd=c("200", "1100"), 
                          name=c("peak1", "peak2")))
    test.rangedData = BED2RangedData(test.bed)
}
\keyword{ misc }
