\name{peakPermTest}
\alias{peakPermTest}
\title{
    Permutation Test for two given peak lists
}
\description{
    Performs a permutation test to seee if there is an association between 
    two given peak lists.
}
\usage{
peakPermTest(peaks1, peaks2, ntimes=100, 
             seed=as.integer(Sys.time()),
             mc.cores=getOption("mc.cores", 2L),
             maxgap=-1L, pool,
             TxDb, bindingDistribution,
             bindingType=c("TSS", "geneEnd"), 
             featureType=c("transcript", "exon"),
             seqn=NA, ...)
}
\arguments{
  \item{peaks1, peaks2}{
  an object of \link[GenomicRanges:GRanges-class]{GRanges}
}
  \item{ntimes}{
  number of permutations
}
  \item{seed}{
  random seed
}
  \item{mc.cores}{
  The number of cores to use. see \link[parallel:mclapply]{mclapply}
}
  \item{maxgap}{
  See \link[IRanges:findOverlaps-methods]{findOverlaps} in the IRanges package for a description 
  of these arguments.
}
  \item{pool}{
  an object of \link{permPool}
}
  \item{TxDb}{
  an object of \link[GenomicFeatures:TxDb-class]{TxDb}
}
  \item{bindingDistribution}{
  an object of \link{bindist}
}
  \item{bindingType}{
  where the peaks should bind, TSS or geneEnd
}
  \item{featureType}{
  what annotation type should be used for detecting the binding distribution.
}
  \item{seqn}{
  default is NA, which means not filter the universe pool for sampling. 
  Otherwise the universe pool will be filtered by the seqnames in seqn.
}
  \item{\dots}{
  further arguments to be passed to \link[regioneR]{numOverlaps}. 
}
}

\value{
A list of class permTestResults. See \link[regioneR]{permTest}
}
\seealso{
    \link{preparePool}, \link{bindist}
}
\references{
Davison, A. C. and Hinkley, D. V. (1997) Bootstrap methods and their 
application, Cambridge University Press, United Kingdom, 156-160
}
\author{
Jianhong Ou
}
\examples{
    path <- system.file("extdata", package="ChIPpeakAnno")
    #files <- dir(path, pattern="[12]_WS170.bed", full.names=TRUE)
    #peaks1 <- toGRanges(files[1], skip=5)
    #peaks2 <- toGRanges(files[2], skip=5)
    #peakPermTest(peaks1, peaks2, TxDb=TxDb.Celegans.UCSC.ce6.ensGene)
    if(interactive()){
        peaks1 <- toGRanges(file.path(path, "MACS2_peaks.xls"), 
                            format="MACS2")
        peaks2 <- toGRanges(file.path(path, "peaks.narrowPeak"), 
                            format="narrowPeak")
        library(TxDb.Hsapiens.UCSC.hg19.knownGene)
        peakPermTest(peaks1, peaks2, 
               TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene, min.pctA=10)
    }
}
\keyword{ misc }
