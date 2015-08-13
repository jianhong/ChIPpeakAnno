\name{permPool-class}
\docType{class}
\alias{permPool}
\alias{permPool-class}
\alias{permPool-method}
\alias{$,permPool-method}
\alias{$<-,permPool-method}

\title{Class \code{"permPool"}}
\description{
  An object of class \code{"permPool"} represents the possible locations to do
  permutation test.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form 
  \code{new("permPool", grs="GRangesList", N="integer")}.
}
\section{Slots}{
  \describe{
    \item{\code{grs}}{object of \code{"GRangesList"} The list of binding ranges}
    \item{\code{N}}{vector of \code{"integer"}, permutation number for each ranges}
  }
}
\section{Methods}{
 \describe{
    \item{$, $<-}{Get or set the slot of \code{\link{permPool}}}
 } 
}
\seealso{
   \link{preparePool}, \link{peakPermTest}
}
\keyword{classes}
