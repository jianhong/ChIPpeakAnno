\name{annotatedPeak}
\alias{annotatedPeak}
\docType{data}
\title{ Annotated Peaks}
\description{
  TSS annotated putative STAT1-binding regions that are identified in 
  un-stimulated cells using ChIP-seq technology (Robertson et al., 2007) 
}
\usage{data(annotatedPeak)}
\format{
   GRanges with slot start holding the start position of the peak, slot end 
   holding the end position of the peak, slot names holding the id of the peak,
   slot strand holding the strands and slot space holding the chromosome 
   location where the peak is located. In addition, the following variables 
   are included.
 \describe{
    \item{\code{feature}}{id of the feature such as ensembl gene ID}
    \item{\code{insideFeature}}{upstream: peak resides upstream of the feature;
    downstream: peak resides downstream of the feature; inside: peak resides 
    inside the feature; overlapStart: peak overlaps with the start of the
    feature; overlapEnd: peak overlaps with the end of the feature; 
    includeFeature: peak include the feature entirely}
    \item{\code{distancetoFeature}}{distance to the nearest feature such as 
    transcription start site}
    \item{\code{start_position}}{start position of the feature such as gene}
    \item{\code{end_position}}{end position of the feature such as the gene}
}
}
\details{
  obtained by 
    data(TSS.human.GRCh37)
  
    data(myPeakList)
  
    annotatePeakInBatch(myPeakList, AnnotationData = TSS.human.GRCh37, 
                        output="b", multiple=F)
}
\examples{
data(annotatedPeak)
str(annotatedPeak) 
if (Sys.getenv("USER")=="jianhongou") {
y = annotatedPeak$distancetoFeature[!is.na(annotatedPeak$distancetoFeature)]
hist(as.numeric(as.character(y)), 
     xlab="Distance To Nearest TSS", main="", breaks=1000, 
ylim=c(0, 50), xlim=c(min(as.numeric(as.character(y)))-100, 
max(as.numeric(as.character(y)))+100))
}
}
\keyword{datasets}
