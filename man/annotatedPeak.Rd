\name{annotatedPeak}
\alias{annotatedPeak}
\docType{data}
\title{ Annotated Peaks}
\description{
  TSS annotated putative STAT1-binding regions that are identified in un-stimulated cells using ChIP-seq technology (Robertson et al., 2007) 
}
\usage{data(annotatedPeak)}
\format{
   RangedData with slot start holding the start position of the peak, slot end holding the end position of the peak, slot rownames holding the id of the peak and slot space holding the chromosome location where the peak is located. In addition, the following variables are included.
 \describe{
    \item{\code{feature}}{id of the feature such as ensembl gene ID}
    \item{\code{insideFeature}}{TRUE for inside the feature and FALSE for between features}
    \item{\code{distancetoFeature}}{distance to the nearest feature such as transcription start site}
    \item{\code{start_position}}{start position of the feature such as gene}
    \item{\code{end_position}}{end position of the feature such as the gene}
    \item{\code{strand}}{1 for positive strand and -1 for negative strand where the feature is located}
}
}
\details{
  obtained by 
	data(TSS.human.NCBI36)
	data(myPeakList)
	annotatePeakInBatch(myPeakList, AnnotationData = TSS.human.NCBI36)
}
\source{
  
}
\references{
  
}
\examples{
data(annotatedPeak)
str(annotatedPeak) 
if (interactive()) {
y = annotatedPeak$distancetoFeature[!is.na(annotatedPeak$distancetoFeature)]
hist(as.numeric(as.character(y)), xlab="Distance To Nearest TSS", main="", breaks=1000, ylim=c(0, 50), xlim=c(min(as.numeric(as.character(y)))-100, max(as.numeric(as.character(y)))+100))
}
}
\keyword{datasets}
