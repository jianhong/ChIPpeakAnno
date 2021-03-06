% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/downstreams.R
\name{downstreams}
\alias{downstreams}
\title{Get downstream coordinates}
\usage{
downstreams(gr, upstream, downstream)
}
\arguments{
\item{gr}{A GenomicRanges object}

\item{upstream, downstream}{non-negative interges.}
}
\value{
A GenomicRanges object
}
\description{
Returns an object of the same type and length as x containing downstream 
ranges. The output range is defined as
}
\details{
(end(x) - upstream) to (end(x) + downstream -1) 

for ranges on the + and * 
strand, and as 

(start(x) - downstream + 1) to (start(x) + downstream)

for ranges on the - strand.

Note that the returned object might contain out-of-bound ranges.
}
\examples{
gr <- GRanges("chr1", IRanges(rep(10, 3), width=6), c("+", "-", "*"))
downstreams(gr, 2, 2)

}
