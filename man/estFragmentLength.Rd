\name{estFragmentLength}
\alias{estFragmentLength}
\title{
    estimate the fragment length
}
\description{
    estimate the fragment length for bam files
}
\usage{
estFragmentLength(bamfiles, index = bamfiles, plot = TRUE, lag.max = 1000, ...)
}
\arguments{
  \item{bamfiles}{
    The file names of the 'BAM' ('SAM' for asBam) files to be processed.
}
  \item{index}{
    The names of the index file of the 'BAM' file being processed; 
    this is given without the '.bai' extension.
}
  \item{plot}{
    logical. If TRUE (the default) the acf is plotted.
}
  \item{lag.max}{
    maximum lag at which to calculate the acf. See \code{\link[stats]{acf}}
}
  \item{\dots}{
    Not used.
}
}
\value{
    numberic vector
}
\author{
    Jianhong Ou
}
\examples{
    if(interactive()){
        library(MMDiffBamSubset)
        WT.AB2 <- WT.AB2()
        Null.AB2 <- Null.AB2()
        Resc.AB2 <- Resc.AB2()
        estFragmentLength(c(WT.AB2, Null.AB2, Resc.AB2))
    }
    
}
\keyword{ misc }
