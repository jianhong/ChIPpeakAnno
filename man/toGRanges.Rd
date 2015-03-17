\name{toGRanges}
\alias{toGRanges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Convert dataset to GRanges
}
\description{
  Convert BED, GFF, RangeData or any user defined dataset to GRanges
}
\usage{
toGRanges(data, format=c("BED", "GFF", "RangedData", "MACS", "others"), 
                      header=FALSE, comment.char="#", colNames=NULL, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  BED, GFF, RangedData or any user defined dataset or their file path.
}
\item{format}{
  data format.  If the data format is set to BED or GFF, please refer to http://genome.ucsc.edu/FAQ/FAQformat#format1 for column order. or MACS output file.
}
\item{header}{
  a logical value indicating whether the file contains the names of the variables as its first line. If missing, the value is determined from the file format: header is set to TRUE if and only if the first row contains one fewer field than the number of columns.
}
\item{comment.char}{
  character: a character vector of length one containing a single character or an empty string. Use "" to turn off the interpretation of comments altogether.
}
  \item{colNames}{
  If the data format is set to "others", colname must be defined. And the colname must contain
  space, start and end. If your columne is names as seqname or chrom, and so on, please rename it as space.
}
  \item{\dots}{
  parameters passed to \link{read.table}
}
}
\value{
 An object of \link[GenomicRanges:GRanges]{GRanges}
}
\author{
  Jianhong Ou
}
\examples{
  rd <- RangedData(IRanges(start = c(967654, 2010897, 2496704),
   end = c(967754, 2010997, 2496804), names = c("Site1", "Site2", "Site3")),
	 space = c("1", "2", "3"), strand=as.integer(1),feature=c("a","b","f"))
  toGRanges(rd, format="RangedData")
}
\keyword{ misc }