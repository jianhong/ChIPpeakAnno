\name{toGRanges}
\alias{toGRanges}
\title{
  Convert dataset to GRanges
}
\description{
  Convert UCSC BED format and its variants, such as GFF, or any user defined dataset such as RangedDate or MACS output file to GRanges
}
\usage{
toGRanges(data, format=c("BED", "GFF", "RangedData", 
                         "MACS", "MACS2", 
                         "narrowPeak", "broadPeak",
                         "others"), 
          header=FALSE, comment.char="#", colNames=NULL, \dots)
}
\arguments{
  \item{data}{
  BED, GFF, RangedData or any user defined dataset or their file path.
}
\item{format}{
  data format.  If the data format is set to BED, GFF, narrowPeak or broadPeak, please refer to http://genome.ucsc.edu/FAQ/FAQformat#format1 for column order. "MACS" is for converting the excel output file from MACS1. "MACS2" is for converting the output file from MACS2.
}
\item{header}{
  A logical value indicating whether the file contains the names of the variables as its first line. If missing, the value is determined from the file format: header is set to TRUE if and only if the first row contains one fewer field than the number of columns.
}
\item{comment.char}{
  character: a character vector of length one containing a single character or an empty string. Use "" to turn off the interpretation of comments altogether.
}
  \item{colNames}{
  If the data format is set to "others", colname must be defined. And the colname must contain
  space, start and end. The column name for the chromosome # should be named as space.
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
  macs <- system.file("extdata", "MACS_peaks.xls", package="ChIPpeakAnno")
  macsOutput <- toGRanges(macs, format="MACS")
}
\keyword{ misc }