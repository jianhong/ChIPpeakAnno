\name{annotatePeakInBatch}
\alias{annotatePeakInBatch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ obtain the distance to the nearest TSS, miRNA, exon et al for a list of peak intervals}
\description{
  obtain the distance to the nearest TSS, miRNA, exon et al for a list of peak locations leveraging IRanges and biomaRt package
}
\usage{
annotatePeakInBatch(myPeakList, mart, featureType = c("TSS", "miRNA","Exon"), 
AnnotationData,output=c("nearestLocation", "overlapping", "both", 
                         "shortestDistance", "inside",
                         "upstream&inside", "inside&downstream",
                         "upstream", "downstream", 
                         "upstreamORdownstream"),
multiple=c(TRUE,FALSE), 
maxgap=0L, PeakLocForDistance = c("start", "middle", "end"), 
FeatureLocForDistance = c("TSS", "middle","start", "end","geneEnd"), 
select=c("all", "first","last","arbitrary"),
ignore.strand=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{myPeakList}{An object of \link[GenomicRanges:GRanges]{GRanges} or \link[IRanges:RangedData]{RangedData}: See example below}
 \item{mart}{used if AnnotationData not supplied, a mart object, see useMart of bioMaRt package for details }
  \item{featureType}{ used if AnnotationData not supplied, TSS, miRNA or exon }
 \item{AnnotationData}{ annotation data obtained from getAnnotation or customized annotation of class RangedData or GRanges containing additional variable:  strand (1 or + for plus strand and -1 or - for minus strand). For example, data(TSS.human.NCBI36),data(TSS.mouse.NCBIM37), data(TSS.rat.RGSC3.4) and data(TSS.zebrafish.Zv8) . If not supplied, then annotation will be obtained from biomaRt automatically using the parameters of mart and featureType}
\item{output}{nearestLocation (default): will output the nearest features calculated as PeakLocForDistance - FeatureLocForDistance;
	overlapping: will output overlapping features with maximum gap specified as maxgap between peak range and feature range;
    shortestDistance: will output nearest features;
	both: will output all the nearest features, in addition, will output any features that overlap the peak that is not the nearest features.
    upstream&inside: will output all upstream and overlapping features with maximum gap.
    inside&downstream: will output all downstream and overlapping features with maximum gap.
    upstream: will output all upstream features with maximum gap.
    downstream: will output all downstream features with maximum gap.
    upstreamORdownstream: will output all upstream features with maximum gap or downstream with maximum gap.
}
\item{multiple}{not applicable when output is nearest. TRUE: output multiple overlapping features for each peak. FALSE: output at most one overlapping feature for each peak. This parameter is kept for backward compatibility, please use select.}
\item{maxgap}{Non-negative integer. Intervals with a separation of maxgap or less are considered to be overlapping}
\item{PeakLocForDistance}{Specify the location of peak for calculating distance,i.e., middle means using middle of the peak to calculate distance to feature, start means using start of the peak to calculate the distance to feature. To be compatible with previous version, by default using start}
\item{FeatureLocForDistance}{Specify the location of feature for calculating distance,i.e., middle means using middle of the feature to calculate distance of peak to feature, start means using start of the feature to calculate the distance to feature, TSS means using start of feature when feature is on plus strand and using end of feature when feature is on minus strand, geneEnd means using end of feature when feature is on plus strand and using start of feature when feature is on minus strand. To be compatible with previous version, by default using TSS}
\item{select}{
all may return multiple overlapping peaks, first will return the first overlapping peak, last will return the last overlapping peak and arbitrary will return one of the overlapping peaks.
}
\item{ignore.strand}{When set to TRUE, the strand information is ignored in the annotation.}
}
\value{
 An object of \link[GenomicRanges:GRanges]{GRanges} or \link[IRanges:RangedData]{RangedData} (depend on what you input) with slot start holding the start position of the peak, slot end holding the end position of the peak, slot space holding the chromosome location where the peak is located, slot rownames holding the id of the peak. In addition, the following variables are included.
    \item{\code{feature}}{id of the feature such as ensembl gene ID}
    \item{\code{insideFeature}}{upstream: peak resides upstream of the feature; downstream: peak resides downstream of the feature; inside: peak resides inside the feature; overlapStart: peak overlaps with the start of the feature; overlapEnd: peak overlaps with the end of the feature; includeFeature: peak include the feature entirely}
    \item{\code{distancetoFeature}}{distance to the nearest feature such as transcription start site.  By default, the distance is calculated as the distance between the start of the binding site and the TSS that is the gene start for genes located on the forward strand and the gene end for genes located on the reverse strand. The user can specify the location of peak and location of feature for calculating this}
    \item{\code{start_position}}{start position of the feature such as gene}
    \item{\code{end_position}}{end position of the feature such as the gene}
    \item{\code{strand}}{1 or + for positive strand and -1 or - for negative strand where the feature is located}
	\item{\code{shortestDistance}}{The shortest distance from either end of peak to either end the feature.
	}
\item{\code{fromOverlappingOrNearest}}{nearest: indicates this feature's start (feature's end for features at minus strand) is closest to the peak start; Overlapping: indicates this feature overlaps with this peak although it is not the nearest feature start
	}
}
\references{ Zhu L.J. et al. (2010) ChIPpeakAnno: a Bioconductor package to annotate ChIP-seq and ChIP-chip data. BMC Bioinformatics 2010, 11:237doi:10.1186/1471-2105-11-237}
\author{ Lihua Julie Zhu, Jianhong Ou}
\seealso{ \link{findOverlappingPeaks}, \link{makeVennDiagram}, \link{addGeneIDs}, \link{peaksNearBDP}, \link{summarizePatternInPeaks}}
\examples{

#if (interactive()){
    ## example 1: annotate myPeakList (RangedData) 
    ## with TSS.human.NCBI36 (RangedData)
    data(myPeakList)
    data(TSS.human.NCBI36)
    annotatedPeak = annotatePeakInBatch(myPeakList[1:6,], 
                                        AnnotationData=TSS.human.NCBI36)
    annotatedPeak
    
    ## example 2: you have a list of transcription factor biding sites from 
    ## literature and are interested in determining the extent of the overlap 
    ## to the list of peaks from your experiment. Prior calling the function 
    ## annotatePeakInBatch, need to represent both dataset as RangedData 
    ## where start is the start of the binding site, end is the end of the 
    ## binding site, names is the name of the binding site, space and strand 
    ## are the chromosome name and strand where the binding site is located.
    
    myexp <- GRanges(seqnames=c(6,6,6,6,5,4,4), 
                     IRanges(start=c(1543200,1557200,1563000,1569800,
                                     167889600,100,1000),
                             end=c(1555199,1560599,1565199,1573799,
                                   167893599,200,1200),
                             names=c("p1","p2","p3","p4","p5","p6", "p7")), 
                     strand="+")
    literature <- GRanges(seqnames=c(6,6,6,6,5,4,4), 
                          IRanges(start=c(1549800,1554400,1565000,1569400,
                                          167888600,120,800),
                                  end=c(1550599,1560799,1565399,1571199,
                                        167888999,140,1400),
                                  names=c("f1","f2","f3","f4","f5","f6","f7")),
                          strand=rep(c("+", "-"), c(5, 2)))
    annotatedPeak1 <- annotatePeakInBatch(myexp, 
                                          AnnotationData=literature)
    pie(table(annotatedPeak1$insideFeature))
    annotatedPeak1
    ### use toGRanges or rtracklayer::import to convert BED or GFF format
    ###  to GRanges before calling annotatePeakInBatch
    test.bed <- data.frame(space=c("4", "6"), 
                           start=c("100", "1000"),
                           end=c("200", "1100"), 
                           name=c("peak1", "peak2"))
    test.GR = toGRanges(test.bed)
    annotatePeakInBatch(test.GR, AnnotationData = literature)
#}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
