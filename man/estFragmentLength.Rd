% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estFragmentLength.R
\name{estFragmentLength}
\alias{estFragmentLength}
\title{estimate the fragment length}
\usage{
estFragmentLength(
  bamfiles,
  index = bamfiles,
  plot = TRUE,
  lag.max = 1000,
  minFragmentSize = 100,
  ...
)
}
\arguments{
\item{bamfiles}{The file names of the 'BAM' ('SAM' for asBam) files to be
processed.}

\item{index}{The names of the index file of the 'BAM' file being processed;
this is given without the '.bai' extension.}

\item{plot}{logical. If TRUE (the default) the acf is plotted.}

\item{lag.max}{maximum lag at which to calculate the acf. See
\code{\link[stats]{acf}}}

\item{minFragmentSize}{minimal fragment size to avoid the phantom peak.}

\item{\dots}{Not used.}
}
\value{
numberic vector
}
\description{
estimate the fragment length for bam files
}
\examples{

    if(interactive() || Sys.getenv("USER")=="jianhongou"){
        path <- system.file("extdata", "reads", package="MMDiffBamSubset")
        if(file.exists(path)){
            WT.AB2 <- file.path(path, "WT_2.bam")
            Null.AB2 <- file.path(path, "Null_2.bam")
            Resc.AB2 <- file.path(path, "Resc_2.bam")
            estFragmentLength(c(WT.AB2, Null.AB2, Resc.AB2))
        }
    }
    

}
\author{
Jianhong Ou
}
\keyword{misc}
