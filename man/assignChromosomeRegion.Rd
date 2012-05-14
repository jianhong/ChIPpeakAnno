\name{assignChromosomeRegion}
\alias{assignChromosomeRegion}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizing peak distribution over exon, intron, enhancer, proximal promoter, 5 prime UTR and 3 prime UTR
}
\description{
Summarizing peak distribution over exon, intron, enhancer, proximal promoter, 5 prime UTR and 3 prime UTR
}
\usage{
assignChromosomeRegion(peaks.RD, exon, TSS, utr5, utr3, proximal.promoter.cutoff = 1000, immediate.downstream.cutoff = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peaks.RD}{
peaks in RangedData: See example below
}
  \item{exon}{
exon data obtained from getAnnotation or customized annotation of class RangedData containing additional variable:  strand (1 or + for plus strand and -1 or - for minus strand)
}
  \item{TSS}{
TSS data obtained from getAnnotation or customized annotation of class RangedData containing additional variable:  strand (1 or + for plus strand and -1 or - for minus strand). For example, data(TSS.human.NCBI36),data(TSS.mouse.NCBIM37), data(TSS.rat.RGSC3.4) and data(TSS.zebrafish.Zv8). 
}
  \item{utr5}{
5 prime UTR data obtained from getAnnotation or customized annotation of class RangedData containing additional variable:  strand (1 or + for plus strand and -1 or - for minus strand).
}
  \item{utr3}{
3 prime UTR data obtained from getAnnotation or customized annotation of class RangedData containing additional variable:  strand (1 or + for plus strand and -1 or - for minus strand). 
}
  \item{proximal.promoter.cutoff}{
Specify the cutoff in bases to be classified as proximal promoter region. Peaks that reside within proximal.promoter.cutoff upstream from or overlap with transcription start site are classified as proximal promoters. Peaks that reside upstream over proximal.promoter.cutoff from gene start are classified as enhancers. The default is 1000 bases.

}
  \item{immediate.downstream.cutoff}{
Specify the cutoff in bases to be classified as immediate downstream. Peaks that reside within immediate.downstream.cutoff downstream of gene end but not overlap 3’ UTR are classified as immediate downstream. Peaks that reside downstream over immediate.downstreatm.cutoff from gene end are classified as enhancers. The default is 1000 bases.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
 \item{\code{Exon}}{Percent of peaks reside in exon regions.}
    \item{\code{Intron}}{Percent of peaks reside in intron regions.}
    \item{\code{5UTR}}{Percent of peaks reside in 5 prime UTR regions.}
    \item{\code{3UTR}}{Percent of peaks reside in 3 prime UTR regions.}
   \item{\code{Proximal Promoter}}{Percent of peaks reside in proximal promoter regions.}
   \item{\code{Immediate Downstream}}{Percent of peaks reside in immediate downstream regions.}
   \item{\code{Enhancer}}{Percent of peaks reside in enhancer regions.}
}
\references{ Zhu L.J. et al. (2010) ChIPpeakAnno: a Bioconductor package to annotate ChIP-seq and ChIP-chip data. BMC Bioinformatics 2010, 11:237doi:10.1186/1471-2105-11-237}
\author{ Lihua Julie Zhu }
\note{ 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{annotatePeakInBatch, findOverlappingPeaks,getEnrichedGO, getEnrichedPATH, makeVennDiagram,addGeneIDs, peaksNearBDP,summarizePatternInPeaks}

\examples{
if (interactive())
{
	library(ChIPpeakAnno)
	data(myPeakList)
	mart<-useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")
	TSS <- getAnnotation(mart, featureType ="TSS")
	utr5 = getAnnotation(mart, featureType = "5utr")
	utr3 = getAnnotation(mart, featureType="3utr")
	exon = getAnnotation(mart, featureType="Exon")
	Feature.distribution = assignChromosomeRegion(myPeakList, exon, TSS, utr5, utr3)
	barplot(unlist(Feature.distribution),cex.names=1,  xlab="Chromosome Region", ylab="Percent Binding Sites", main="")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
