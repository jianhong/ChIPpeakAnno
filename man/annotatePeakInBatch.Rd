\name{annotatePeakInBatch}
\alias{annotatePeakInBatch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ obtain the distance to the nearest TSS, miRNA, exon et al for a list of peak intervals}
\description{
  obtain the distance to the nearest TSS, miRNA, exon et al for a list of peak locations leveraging IRanges and biomaRt package
}
\usage{
annotatePeakInBatch(myPeakList, mart, featureType = c("TSS", "miRNA","Exon"), AnnotationData)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{myPeakList}{RangedData: See example below}
  \item{mart}{used if AnnotationData not supplied, a mart object, see useMart of bioMaRt package for details }
   \item{featureType}{ used if AnnotationData not supplied, TSS, miRNA or exon }
  \item{AnnotationData}{ annotation data obtained from getAnnotation or customized annotation of class RangedData containing additional variable:  strand (1 or + for plus strand and -1 or - for minus strand). For example, data(TSS.human.NCBI36),data(TSS.mouse.NCBIM37), data(GO.rat.RGSC3.4) and data(TSS.zebrafish.Zv8) . If not supplied, then annotation will be obtained from biomaRt automatically using the parameters of mart and featureType}
}
\details{
  }
\value{
 RangedData with slot start holding the start position of the peak, slot end holding the end position of the peak, slot space holding the chromosome location where the peak is located, slot rownames holding the id of the peak. In addition, the following variables are included.
    \item{\code{feature}}{id of the feature such as ensembl gene ID}
    \item{\code{insideFeature}}{upstream: peak resides upstream of the feature; downstream: peak resides downstream of the feature; inside: peak resides inside the feature; overlapStart: peak overlaps with the start of the feature; overlapEnd: peak overlaps with the end of the feature; includeFeature: peak include the feature entirely}
    \item{\code{distancetoFeature}}{distance to the nearest feature such as transcription start site. The distance is calculated as the distance between the start of the binding site and the TSS that is the gene start for genes located on the forward strand and the gene end for genes located on the reverse strand}
    \item{\code{start_position}}{start position of the feature such as gene}
    \item{\code{end_position}}{end position of the feature such as the gene}
    \item{\code{strand}}{1 or + for positive strand and -1 or - for negative strand where the feature is located}
}
\references{ Durinck S. et al. (2005) BioMart and Bioconductor: a powerful link between biological biomarts and microarray data analysis. Bioinformatics, 21, 3439-3440.  }
\author{ Lihua Julie Zhu }
\note{ 
}
\seealso{ findOverlappingPeaks, makeVennDiagram }
\examples{

if (interactive())
{
	## example 1: annotate myPeakList (RangedData) with TSS.human.NCBI36 (RangedData)
	data(myPeakList)
	data(TSS.human.NCBI36)
	annotatedPeak = annotatePeakInBatch(myPeakList[1:6,], AnnotationData=TSS.human.NCBI36)
	as.data.frame(annotatedPeak)
	## example 2: you have a list of transcription factor biding sites from literature and are interested in 
	## determining the extent of the overlap to the list of peaks from your experiment
	## Prior calling the function annotatePeakInBatch, need to represent both dataset as RangedData where start is the start 
	## of the binding site, end is the end of the binding site, names is the name of the binding site, 
	## space and strand are the chromosome name and strand where the binding site is located.

	myexp =  RangedData(IRanges(start=c(1543200,1557200,1563000,1569800,167889600,100,1000),end=c(1555199,1560599,1565199,1573799,167893599,200,1200),names=c("p1","p2","p3","p4","p5","p6", "p7")),strand=as.integer(1),space=c(6,6,6,6,5,4,4))	literature = RangedData(IRanges(start=c(1549800,1554400,1565000,1569400,167888600,120,800),end=c(1550599,1560799,1565399,1571199,167888999,140,1400),names=c("f1","f2","f3","f4","f5","f6","f7")),strand=c(1,1,1,1,1,-1,-1),space=c(6,6,6,6,5,4,4))	annotatedPeak1= annotatePeakInBatch(myexp, AnnotationData = literature)
	pie(table(as.data.frame(annotatedPeak1)$insideFeature))
	as.data.frame(annotatedPeak1)
	### use BED2RangedData or GFF2RangedData to convert BED format or GFF format to RangedData before calling annotatePeakInBatch
	test.bed = data.frame(cbind(chrom = c("4", "6"), chromStart=c("100", "1000"),chromEnd=c("200", "1100"), name=c("peak1", "peak2")))
	test.rangedData = BED2RangedData(test.bed)
	annotatePeakInBatch(test.rangedData, AnnotationData = literature)
	test.GFF = data.frame(cbind(seqname  = c("chr4", "chr4"), source=rep("Macs", 2), feature=rep("peak", 2), start=c("100", "1000"), end=c("200", "1100"), score=c(60, 26), strand=c(1, 1), frame=c(".", 2), group=c("peak1", "peak2")))
	test.rangedData = GFF2RangedData(test.GFF)
	as.data.frame(annotatePeakInBatch(test.rangedData, AnnotationData = literature))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }