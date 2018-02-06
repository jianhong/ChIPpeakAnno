\name{peaksNearBDP}
\alias{peaksNearBDP}
\title{ obtain the peaks near bi-directional promoters}
\description{
  Obtain the peaks near bi-directional promoters. Also output percent of 
  peaks near bi-directional promoters.
}
\usage{
peaksNearBDP(myPeakList, AnnotationData, MaxDistance=5000L, ...)
}
\arguments{
 \item{myPeakList}{\link[GenomicRanges:GRanges-class]{GRanges} or 
 \link[IRanges:RangedData-class]{RangedData}: See example below}
 \item{AnnotationData}{ annotation data obtained from getAnnotation or 
 customized annotation of class \link[GenomicRanges:GRanges-class]{GRanges} containing additional variable:  strand (1 or + for plus 
 strand and -1 or - for minus strand). For example, data(TSS.human.NCBI36),
 data(TSS.mouse.NCBIM37), data(TSS.rat.RGSC3.4) and data(TSS.zebrafish.Zv8).}
\item{MaxDistance}{Specify the maximum gap allowed between the peak and 
nearest gene}
\item{...}{Not used}
}
\value{
A list of 4
 \item{\code{peaksWithBDP}}{annotated Peaks containing bi-directional promoters.

 GRangesList with slot start holding the start position of the peak, slot end 
 holding the end position of the peak, slot space holding the chromosome 
 location where the peak is located, slot rownames holding the id of the peak. 
 In addition, the following variables are included.

    feature: id of the feature such as ensembl gene ID

    insideFeature: upstream: peak resides upstream of the feature; downstream: 
    peak resides downstream of the feature; inside: peak resides inside the 
    feature; overlapStart: peak overlaps with the start of the feature; 
    overlapEnd: peak overlaps with the end of the feature; includeFeature: 
    peak include the feature entirely.

    distancetoFeature: distance to the nearest feature such as transcription 
    start site.  By default, the distance is calculated as the distance between
    the start of the binding site and the TSS that is the gene start for genes
    located on the forward strand and the gene end for genes located on the 
    reverse strand. The user can specify the location of peak and location of 
    feature for calculating this

    feature_range: start and end position of the feature such as gene

    feature_strand: 1 or + for positive strand and -1 or - for negative strand where 
    the feature is located
}
\item{\code{percentPeaksWithBDP}}{The percent of input peaks containing 
bi-directional promoters}
\item{\code{n.peaks}}{The total number of input peaks}
\item{\code{n.peaksWithBDP}}{The # of input peaks containing bi-directional 
promoters}
}
\references{ Zhu L.J. et al. (2010) ChIPpeakAnno: a Bioconductor package to 
annotate ChIP-seq and ChIP-chip data. BMC Bioinformatics 2010, 
11:237doi:10.1186/1471-2105-11-237}
\author{ Lihua Julie Zhu, Jianhong Ou}
\seealso{annotatePeakInBatch, findOverlappingPeaks, makeVennDiagram }
\examples{

if (Sys.getenv("USER")=="jianhongou")
{
    data(myPeakList)
    data(TSS.human.NCBI36)
    seqlevelsStyle(TSS.human.NCBI36) <- seqlevelsStyle(myPeakList)
    annotatedBDP = peaksNearBDP(myPeakList[1:6,], 
                                AnnotationData=TSS.human.NCBI36,
                                MaxDistance=5000,
                                PeakLocForDistance =  "middle", 
                                FeatureLocForDistance = "TSS")
    c(annotatedBDP$percentPeaksWithBDP, annotatedBDP$n.peaks, 
      annotatedBDP$n.peaksWithBDP)
}
}
\keyword{ misc }
