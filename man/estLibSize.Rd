\name{estLibSize}
\alias{estLibSize}
\title{
    estimate the library size
}
\description{
    estimate the library size of bam files
}
\usage{
estLibSize(bamfiles, index = bamfiles, ...)
}
\arguments{
  \item{bamfiles}{
    The file names of the 'BAM' ('SAM' for asBam) files to be processed.
}
  \item{index}{
    The names of the index file of the 'BAM' file being processed; 
    this is given without the '.bai' extension.
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
        estLibSize(c(WT.AB2, Null.AB2, Resc.AB2))
    }
}

\keyword{ misc }
