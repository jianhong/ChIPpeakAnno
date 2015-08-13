\name{bindist-class}
\docType{class}
\alias{bindist}
\alias{bindist-class}
\alias{bindist-method}
\alias{$,bindist-method}
\alias{$<-,bindist-method}

\title{Class \code{"bindist"}}
\description{
  An object of class \code{"bindist"} represents the relevant fixed-width range
  of binding site from the feature and number of possible binding site 
  in each range.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form 
  \code{new("bindist", counts="integer",
                            mids="integer",
                            halfBinSize="integer",
                            bindingType="character",
                            featureType="character")}.
}
\section{Slots}{
  \describe{
    \item{\code{counts}}{vector of \code{"integer"} The count number in each binding range}
    \item{\code{mids}}{vector of \code{"integer"} The center of each range relevant to feature}
    \item{\code{halfBinSize}}{\code{"integer"}, length must be 1. the fixed half-width of each binding range}
    \item{\code{bindingType}}{a \code{"character"}. could be "TSS", "geneEnd"}
    \item{\code{featureType}}{a \code{"character"}. could be "transcript", "exon"}
  }
}
\section{Methods}{
 \describe{
    \item{$, $<-}{Get or set the slot of \code{\link{bindist}}}
 } 
}
\seealso{
   \link{preparePool}, \link{peakPermTest}
}
\keyword{classes}
