% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/myPeakList.R
\docType{data}
\name{myPeakList}
\alias{myPeakList}
\title{An example GRanges object representing a ChIP-seq peak dataset}
\format{
GRanges with slot rownames containing the ID of peak as character,
slot start containing the start position of the peak, slot end containing
the end position of the peak and seqnames containing the chromosome where
the peak is located.
}
\source{
Robertson G, Hirst M, Bainbridge M, Bilenky M, Zhao Y, et al. (2007)
Genome-wide profiles of STAT1 DNA association using chromatin
immunoprecipitation and massively parallel sequencing. Nat Methods 4:651-7
}
\usage{
myPeakList
}
\description{
the putative STAT1-binding regions identified in un-stimulated cells using
ChIP-seq technology (Robertson et al., 2007)
}
\examples{

data(myPeakList)
slotNames(myPeakList)
}
\keyword{datasets}
