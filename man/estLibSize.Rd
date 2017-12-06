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
    if(Sys.getenv("USER")=="jianhongou"){
        path <- system.file("extdata", "reads", package="MMDiffBamSubset")
        if(file.exists(path)){
            WT.AB2 <- file.path(path, "WT_2.bam")
            Null.AB2 <- file.path(path, "Null_2.bam")
            Resc.AB2 <- file.path(path, "Resc_2.bam")
            estLibSize(c(WT.AB2, Null.AB2, Resc.AB2))
        }
    }
}

\keyword{ misc }
