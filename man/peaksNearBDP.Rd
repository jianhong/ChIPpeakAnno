\name{peaksNearBDP}
\alias{peaksNearBDP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ obtain the peaks near bi-directional promoters}
\description{
  Obtain the peaks near bi-directional promoters. Also output % of peaks near bi-directional promoters.
}
\usage{
peaksNearBDP(myPeakList, mart,AnnotationData, MaxDistance=5000,PeakLocForDistance = c("start", "middle", "end"), 
FeatureLocForDistance = c("TSS", "middle","start", "end","geneEnd"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{myPeakList}{RangedData: See example below}
 \item{mart}{used if AnnotationData not supplied, a mart object, see useMart of bioMaRt package for details }
 \item{AnnotationData}{ annotation data obtained from getAnnotation or customized annotation of class RangedData containing additional variable:  strand (1 or + for plus strand and -1 or - for minus strand). For example, data(TSS.human.NCBI36),data(TSS.mouse.NCBIM37), data(TSS.rat.RGSC3.4) and data(TSS.zebrafish.Zv8) . If not supplied, then annotation will be obtained from biomaRt automatically using the parameters of mart and featureType TSS}
\item{MaxDistance}{Specify the maximum gap allowed between the peak and nearest gene}
\item{PeakLocForDistance}{Specify the location of peak for calculating distance,i.e., middle means using middle of the peak to calculate distance to feature, start means using start of the peak to calculate the distance to feature. To be compatible with previous version, by default using start}
\item{FeatureLocForDistance}{Specify the location of feature for calculating distance,i.e., middle means using middle of the feature to calculate distance of peak to feature, start means using start of the feature to calculate the distance to feature, TSS means using start of feature when feature is on plus strand and using end of feature when feature is on minus strand, geneEnd means using end of feature when feature is on plus strand and using start of feature when feature is on minus strand. To be compatible with previous version, by default using TSS}
}
\details{
  }
\value{
A list of 4
 \item{\code{peaksWithBDP}}{annotated Peaks containing bi-directional promoters.

 RangedData with slot start holding the start position of the peak, slot end holding the end position of the peak, slot space holding the chromosome location where the peak is located, slot rownames holding the id of the peak. In addition, the following variables are included.

    feature: id of the feature such as ensembl gene ID

    insideFeature: upstream: peak resides upstream of the feature; downstream: peak resides downstream of the feature; inside: peak resides inside the feature; overlapStart: peak overlaps with the start of the feature; overlapEnd: peak overlaps with the end of the feature; includeFeature: peak include the feature entirely.

    distancetoFeature: distance to the nearest feature such as transcription start site.  By default, the distance is calculated as the distance between the start of the binding site and the TSS that is the gene start for genes located on the forward strand and the gene end for genes located on the reverse strand. The user can specify the location of peak and location of feature for calculating this

    start_position: start position of the feature such as gene

    end_position: end position of the feature such as the gene

    strand: 1 or + for positive strand and -1 or - for negative strand where the feature is located

    shortestDistance: The shortest distance from either end of peak to either end the feature

    fromOverlappingOrNearest: NearestStart: indicates this PeakLocForDistance is closest to the 
FeatureLocForDistance
}
\item{\code{percentPeaksWithBDP}}{The percent of input peaks containing bi-directional promoters}
\item{\code{n.peaks}}{The total number of input peaks}
\item{\code{n.peaksWithBDP}}{The # of input peaks containing bi-directional promoters}
}
\references{ Zhu L.J. et al. (2010) ChIPpeakAnno: a Bioconductor package to annotate ChIP-seq and ChIP-chip data. BMC Bioinformatics 2010, 11:237doi:10.1186/1471-2105-11-237}
\author{ Lihua Julie Zhu }
\note{ 
}
\seealso{annotatePeakInBatch, findOverlappingPeaks, makeVennDiagram }
\examples{

if (interactive())
{
	data(myPeakList)
	data(TSS.human.NCBI36)
	annotatedBDP = peaksNearBDP(myPeakList[1:6,], AnnotationData=TSS.human.NCBI36,
	MaxDistance=5000,PeakLocForDistance =  "middle", 
	FeatureLocForDistance = "TSS")
	c(annotatedBDP$percentPeaksWithBDP, annotatedBDP$n.peaks, annotatedBDP$n.peaksWithBDP)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
