% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binOverRegions.R
\name{plotBinOverRegions}
\alias{plotBinOverRegions}
\title{plot the coverage of regions}
\usage{
plotBinOverRegions(dat, ...)
}
\arguments{
\item{dat}{A list of matrix which indicate the coverage of regions per bin}

\item{...}{Parameters could be used by \link[graphics]{matplot}}
}
\description{
plot the output of \link{binOverRegions} or \link{binOverGene}
}
\examples{

if(interactive()){
path <- system.file("extdata", package="ChIPpeakAnno")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(rtracklayer)
files <- dir(path, "bigWig")
if(.Platform$OS.type != "windows"){
cvglists <- lapply(file.path(path, files), import,
                   format="BigWig", as="RleList")
names(cvglists) <- sub(".bigWig", "", files)
d <- binOverGene(cvglists, TxDb.Hsapiens.UCSC.hg19.knownGene)
plotBinOverRegions(d)
}
}

}
\seealso{
\link{binOverRegions}, \link{binOverGene}
}
\author{
Jianhong Ou
}
