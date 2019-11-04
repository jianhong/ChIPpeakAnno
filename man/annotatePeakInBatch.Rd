\name{annotatePeakInBatch}
\alias{annotatePeakInBatch}
\title{ Obtain the distance to the nearest TSS, miRNA, and/or exon for a list 
of peaks}
\description{
  Obtain the distance to the nearest TSS, miRNA, exon et al for a list of peak 
  locations leveraging IRanges and biomaRt package
}
\usage{
annotatePeakInBatch(myPeakList, mart, featureType = c("TSS", "miRNA","Exon"), 
AnnotationData, output=c("nearestLocation", "overlapping", "both", 
                         "shortestDistance", "inside",
                         "upstream&inside", "inside&downstream",
                         "upstream", "downstream", 
                         "upstreamORdownstream",
                         "nearestBiDirectionalPromoters"),
multiple=c(TRUE,FALSE), 
maxgap=-1L, PeakLocForDistance=c("start", "middle", "end"),
FeatureLocForDistance=c("TSS", "middle","start", "end","geneEnd"), 
select=c("all", "first","last","arbitrary"),
ignore.strand=TRUE, bindingRegion=NULL, ...)
}
\arguments{
 \item{myPeakList}{A \link[GenomicRanges:GRanges-class]{GRanges} object }
 \item{mart}{A mart object, used if AnnotationData is not supplied, see useMart 
 of bioMaRt package for details }
 \item{featureType}{ A charcter vector used with mart argument if 
 AnnotationData is not supplied; it's value is "TSS"", "miRNA"" or "Exon" }
 \item{AnnotationData}{ A \link[GenomicRanges:GRanges-class]{GRanges} or \link{annoGR} oject. 
 It can be obtained from function getAnnotation or customized annotation of 
 class GRanges containing additional variable:  strand (1 or + for plus strand 
 and -1 or - for minus strand). Pre-compliled annotations, such as 
 TSS.human.NCBI36, TSS.mouse.NCBIM37, TSS.rat.RGSC3.4 and TSS.zebrafish.Zv8, 
 are provided by this package (attach them with data() function). Another 
 method to provide annotation data is to obtain through biomaRt real time by 
 using the parameters of mart and featureType}
\item{output}{
    \describe{
        \item{nearestLocation (default)}{will output the nearest features calculated as PeakLoc - FeatureLocForDistance}
        \item{overlapping}{will output overlapping features with maximum gap specified as maxgap between peak range and feature range}
        \item{shortestDistance}{will output nearest features}
        \item{both}{will output all the nearest features, in addition, will output any  features that overlap the peak that is not the nearest features}
        \item{upstream&inside}{will output all upstream and overlapping features with maximum gap}
        \item{inside&downstream}{will output all downstream and overlapping features 
with maximum gap}
        \item{upstream}{will output all upstream features with maximum gap.}
        \item{downstream}{will output all downstream features with maximum gap.}
        \item{upstreamORdownstream}{will output all upstream features with maximum gap or downstream with maximum gap}
        \item{nearestBiDirectionalPromoters}{will use \link{annoPeaks} to annotate peaks. Nearest promoters from both direction of the peaks (strand is considered). It will report bidirectional promoters if there are promoters in both directions in the given region (defined by bindingRegion). Otherwise, it will report the closest promoter in one direction.}
    }
}
\item{multiple}{Not applicable when output is nearest. TRUE: output multiple 
overlapping features for each peak. FALSE: output at most one overlapping 
feature for each peak. This parameter is kept for backward compatibility, 
please use select.}
\item{maxgap}{The maximum \emph{gap} that is allowed between 2 ranges for the
ranges to be considered as overlapping. The \emph{gap} between 2 ranges
is the number of positions that separate them. The \emph{gap} between
2 adjacent ranges is 0. By convention when one range has its start or
end strictly inside the other (i.e. non-disjoint ranges), the \emph{gap}
is considered to be -1.}
\item{PeakLocForDistance}{Specify the location of peak for calculating 
distance,i.e., middle means using middle of the peak to calculate distance 
to feature, start means using start of the peak to calculate the distance 
to feature. To be compatible with previous version, by default using start}
\item{FeatureLocForDistance}{Specify the location of feature for calculating 
distance,i.e., middle means using middle of the feature to calculate distance 
of peak to feature, start means using start of the feature to calculate the
distance to feature, TSS means using start of feature when feature is on plus 
strand and using end of feature when feature is on minus strand, geneEnd means 
using end of feature when feature is on plus strand and using start of feature
when feature is on minus strand. To be compatible with previous version, by 
default using TSS}
\item{select}{
"all" may return multiple overlapping peaks, "first" will return the first 
overlapping peak, "last" will return the last overlapping peak and "arbitrary"
will return one of the overlapping peaks.
}
\item{ignore.strand}{When set to TRUE, the strand information is ignored in 
the annotation.}
\item{bindingRegion}{Annotation range used for \link{annoPeaks}, which is a vector with two integer values, default to c (-5000, 5000). The first one must be no bigger than 0. And the sec ond one must be no less than 1. Once bindingRegion is defined,
annotation will based on \link{annoPeaks}. Here is how to use it together with the parameter output and FeatureLocForDistance.
    \itemize{
        \item To obtain peaks with nearest bi-directional promoters within 5kb upstream and 3kb downstream of TSS, set  output = "nearestBiDirectionalPromoters" and bindingRegion = c(-5000, 3000)
        \item To obtain peaks within 5kb upstream and up to 3kb downstream of TSS within the gene body, set output="overlapping", FeatureLocForDistance="TSS" and bindingRegion = c(-5000, 3000)
        \item To obtain  peaks  up to 5kb upstream within the gene body and 3kb downstream of gene/Exon End, set output="overlapping", FeatureLocForDistance="geneEnd" and bindingRegion = c(-5000, 3000)
        \item To obtain peaks from 5kb upstream to 3kb downstream of genes/Exons, set output="overlapping",  bindingType = "fullRange" and bindingRegion = c(-5000, 3000)
    }
For details, see \link{annoPeaks}.}
\item{...}{Parameters could be passed to \link{annoPeaks}}
}
\value{
An object of \link[GenomicRanges:GRanges-class]{GRanges} with slot start holding 
the start position of the peak, slot end holding the end position of the peak, 
slot space holding the chromosome location where the peak is located, slot 
rownames holding the id of the peak. In addition, the following variables 
are included.
    \item{\code{feature}}{id of the feature such as ensembl gene ID}
    \item{\code{insideFeature}}{upstream: peak resides upstream of the 
    feature; downstream: peak resides downstream of the feature; inside: 
    peak resides inside the feature; overlapStart: peak overlaps with the 
    start of the feature; overlapEnd: peak overlaps with the end of the 
    feature; includeFeature: peak include the feature entirely}
    \item{\code{distancetoFeature}}{distance to the nearest feature such as 
    transcription start site.  By default, the distance is calculated as the
    distance between the start of the binding site and the TSS that is the gene
    start for genes located on the forward strand and the gene end for genes 
    located on the reverse strand. The user can specify the location of peak 
    and location of feature for calculating this}
    \item{\code{start_position}}{start position of the feature such as gene}
    \item{\code{end_position}}{end position of the feature such as the gene}
    \item{\code{strand}}{1 or + for positive strand and -1 or - for negative 
    strand where the feature is located}
    \item{\code{shortestDistance}}{The shortest distance from either end of 
    peak to either end the feature.
    }
\item{\code{fromOverlappingOrNearest}}{nearest: indicates this feature's start 
(feature's end for features at minus strand) is closest to the peak start;
Overlapping: indicates this feature overlaps with this peak although it is
not the nearest feature start
    }
}
\references{ 
1. Zhu L.J. et al. (2010) ChIPpeakAnno: a Bioconductor package to annotate 
ChIP-seq and ChIP-chip data. BMC Bioinformatics 2010, 
11:237doi:10.1186/1471-2105-11-237

2. Zhu L (2013). "Integrative analysis of ChIP-chip and ChIP-seq dataset." 
In Lee T and Luk ACS (eds.), Tilling Arrays, volume 1067, chapter 4, pp. -19. 
Humana Press. http://dx.doi.org/10.1007/978-1-62703-607-8_8
} 

\author{ Lihua Julie Zhu, Jianhong Ou}
\seealso{ \link{getAnnotation}, \link{findOverlappingPeaks}, 
\link{makeVennDiagram}, \link{addGeneIDs}, \link{peaksNearBDP}, 
\link{summarizePatternInPeaks}, \link{annoGR}, \link{annoPeaks}}
\examples{

#if (interactive()){
    ## example 1: annotate myPeakList by TxDb or EnsDb.
    data(myPeakList)
    library(EnsDb.Hsapiens.v75)
    annoData <- annoGR(EnsDb.Hsapiens.v75)
    annotatePeak = annotatePeakInBatch(myPeakList[1:6], AnnotationData=annoData)
    annotatePeak
    
    ## example 2: annotate myPeakList (GRanges) 
    ## with TSS.human.NCBI36 (Granges)
    data(TSS.human.NCBI36)
    annotatedPeak = annotatePeakInBatch(myPeakList[1:6], 
                                        AnnotationData=TSS.human.NCBI36)
    annotatedPeak
    
    ## example 3: you have a list of transcription factor biding sites from 
    ## literature and are interested in determining the extent of the overlap 
    ## to the list of peaks from your experiment. Prior calling the function 
    ## annotatePeakInBatch, need to represent both dataset as GRanges 
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
\keyword{ misc }
