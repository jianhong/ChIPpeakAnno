\name{oligoFrequency}
\alias{oligoFrequency}
\title{
 get the oligonucleotide frequency
}
\description{
 Prepare the oligonucleotide frequency for given Markov order. 
}
\usage{
oligoFrequency(sequence, MarkovOrder = 3L, last = 1e+06)
}
\arguments{
  \item{sequence}{
    The sequences packaged in DNAStringSet, DNAString object or output of 
    function \link{getAllPeakSequence}.
}
  \item{MarkovOrder}{
    Markov order.
}
  \item{last}{
    The sequence size to be analyzed.
}
}

\value{
    A numeric vector.
}
\author{
Jianhong Ou
}
\seealso{
See Also as \code{\link{oligoSummary}}
}
\examples{
    oligoFrequency(DNAString("AATTCGACGTACAGATGACTAGACT"))
}
\keyword{ misc }
