\name{binOverFeature}
\alias{binOverFeature}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
peak aggregation over bins from TSS
}
\description{
peak aggregation over bins from feature sites.
}
\usage{
binOverFeature(..., annotationData=GRanges(),
                    select=c("all", "nearest"),
                    radius=5000L, nbins=50L,
                    minGeneLen=1L, aroundGene=FALSE, mbins=nbins, 
                    featureSite=c("FeatureStart", "FeatureEnd", "bothEnd"),
                    PeakLocForDistance=c("all", "end","start","middle"), 
                    FUN=sum, xlab, ylab, main)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
    objects of GRanges to be analyzed
}
  \item{annotationData}{
    an object of GRanges for annotation
}
  \item{select}{
    annotate the peaks to all features or the nearest one
}
  \item{radius}{
    radius of the longest distance to feature site
}
  \item{nbins}{
    number of bins
}
  \item{minGeneLen}{
    minimal gene length
}
  \item{aroundGene}{
    count peaks around features or a give site of the features
}
  \item{mbins}{
    if aroundGene set as TRUE, the number of bins intra-feature
}
  \item{featureSite}{
    which site of features should be used for distance calculation
}
  \item{PeakLocForDistance}{
    which site of peaks should be used for distance calculation
}
  \item{FUN}{
    the function to be used for score calculation
}
  \item{xlab}{
    titles for each x axis
}
  \item{ylab}{
    titles for each y axis
}
  \item{main}{
    overall titles for each plot
}
}
\value{
 an object of data.frame with bin values.
}
\author{
Jianhong Ou
}
\examples{
bed <- system.file("extdata", "MACS_output.bed", package="ChIPpeakAnno")
gr1 <- toGRanges(bed, format="BED", header=FALSE)
data(TSS.human.GRCh37)
binOverFeature(gr1, annotationData=TSS.human.GRCh37,
               radius=5000, nbins=10, FUN=length)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
