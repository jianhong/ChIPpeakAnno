\name{GFF2RangedData}
\alias{GFF2RangedData}
\alias{GFF2RangedData-deprecated}
\title{ Convert GFF format to RangedData }
\description{
  Convert GFF format to RangedData. This function will be depreciated. 
  Use function toGRanges instead.
 }
\usage{
GFF2RangedData(data.GFF,header=FALSE, ...)
}
\arguments{
  \item{data.GFF}{ GFF format data frame or GFF file name, please refer to 
  http://genome.ucsc.edu/FAQ/FAQformat#format3 for details }
  \item{header}{TRUE or FALSE, default to FALSE, indicates whether data.GFF 
  file has GFF header} 
  \item{...}{any parameter need to be passed into read.delim function}
}
\value{
   RangedData with slot start holding the start position of the feature, slot 
   end holding the end position of the feature, slot names holding the id of 
   the feature, slot space holding the chromosome location where the feature 
   is located. In addition, the following variables are included.
    \item{\code{strand}}{1 for positive strand and -1 for negative strand where
    the feature is located.}
}
\author{ Lihua Julie Zhu }
\note{ 
For converting the peakList in GFF format to RangedData before calling 
annotatePeakInBatch function
}
\examples{
  test.GFF = data.frame(cbind(seqname  = c("chr1", "chr2"), 
  source=rep("Macs", 2), 
  feature=rep("peak", 2), 
  start=c("100", "1000"), 
  end=c("200", "1100"), 
  score=c(60, 26), 
  strand=c(1, -1), 
  frame=c(".", 2), 
  group=c("peak1", "peak2")))
  test.rangedData = GFF2RangedData(test.GFF)
}
\keyword{ misc }
