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
assignChromosomeRegion(peaks.RD, exon, TSS, utr5, utr3, proximal.promoter.cutoff=1000L, immediate.downstream.cutoff=1000L, nucleotideLevel=FALSE, precedence=NULL, TxDb=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peaks.RD}{
peaks in RangedData or GRanges: See example below
}
  \item{exon}{
exon data obtained from getAnnotation or customized annotation of class RangedData containing additional variable:  strand (1 or + for plus strand and -1 or - for minus strand). Will not use anymore! use \code{\link[GenomicFeatures]{TxDb}} instead.
}
  \item{TSS}{
TSS data obtained from getAnnotation or customized annotation of class RangedData containing additional variable:  strand (1 or + for plus strand and -1 or - for minus strand). For example, data(TSS.human.NCBI36),data(TSS.mouse.NCBIM37), data(TSS.rat.RGSC3.4) and data(TSS.zebrafish.Zv8). Will not use anymore! use \code{\link[GenomicFeatures]{TxDb}} instead.
}
  \item{utr5}{
5 prime UTR data obtained from getAnnotation or customized annotation of class RangedData containing additional variable:  strand (1 or + for plus strand and -1 or - for minus strand). Will not use anymore! use \code{\link[GenomicFeatures]{TxDb}} instead.
}
  \item{utr3}{
3 prime UTR data obtained from getAnnotation or customized annotation of class RangedData containing additional variable:  strand (1 or + for plus strand and -1 or - for minus strand). Will not use anymore! use \code{\link[GenomicFeatures]{TxDb}} instead.
}
  \item{proximal.promoter.cutoff}{
Specify the cutoff in bases to be classified as proximal promoter region. Peaks that reside within proximal.promoter.cutoff upstream from or overlap with transcription start site are classified as proximal promoters. Peaks that reside upstream over proximal.promoter.cutoff from gene start are classified as enhancers. The default is 1000 bases.

}
  \item{immediate.downstream.cutoff}{
Specify the cutoff in bases to be classified as immediate downstream. Peaks that reside within immediate.downstream.cutoff downstream of gene end but not overlap 3 prime UTR are classified as immediate downstream. Peaks that reside downstream over immediate.downstreatm.cutoff from gene end are classified as enhancers. The default is 1000 bases.
}
  \item{nucleotideLevel}{NucleotideLevel (TRUE or FALSE) to allow both peak centric and nucleotide centric view. Default=FALSE}
  \item{precedence}{If no precedence specified,
  double count will be enabled, which means that if a peak overlap with both
promoter and 5'UTR, then both promoter and 5'UTR will be incremented. If a
precedence order is specified, for example, if promoter is specified before
5'UTR, then only promoter will be incremented for the same example. Default=NULL}
  \item{TxDb}{an object of \code{\link[GenomicFeatures]{TxDb}}}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
\item{\code{jaccard}}{Jaccard Index}
 \item{\code{Exons}}{Percent of peaks reside in exon regions.}
    \item{\code{Introns}}{Percent of peaks reside in intron regions.}
    \item{\code{fiveUTRs}}{Percent of peaks reside in 5 prime UTR regions.}
    \item{\code{threeUTRs}}{Percent of peaks reside in 3 prime UTR regions.}
   \item{\code{Promoter}}{Percent of peaks reside in proximal promoter regions.}
   \item{\code{ImmediateDownstream}}{Percent of peaks reside in immediate downstream regions.}
   \item{\code{Intergenic.Region}}{Percent of peaks reside in intergenic regions except above regions.}
}
\references{ Zhu L.J. et al. (2010) ChIPpeakAnno: a Bioconductor package to annotate ChIP-seq and ChIP-chip data. BMC Bioinformatics 2010, 11:237doi:10.1186/1471-2105-11-237}
\author{ Jianhong Ou, Lihua Julie Zhu }
\note{ 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{annotatePeakInBatch, findOverlappingPeaks,getEnrichedGO, getEnrichedPATH, makeVennDiagram,addGeneIDs, peaksNearBDP,summarizePatternInPeaks}

\examples{
if (interactive())
{
	library(ChIPpeakAnno)
	##Display the list of genomes available at UCSC:
  #library(rtracklayer)
  #ucscGenomes()[, "db"]
  ## Display the list of Tracks supported by makeTranscriptDbFromUCSC()
  #supportedUCSCtables()
  ##Retrieving a full transcript dataset for Human from UCSC
  TxDb <- makeTranscriptDbFromUCSC(genome="hg19", tablename="ensGene")
  exons <- exons(TxDb, columns=NULL)
  fiveUTRs <- unique(unlist(fiveUTRsByTranscript(TxDb)))
  Feature.distribution <- assignChromosomeRegion(exons, nucleotideLevel=TRUE, TxDb=TxDb)
  barplot(Feature.distribution$percentage)
  assignChromosomeRegion(fiveUTRs, nucleotideLevel=FALSE, TxDb=TxDb)
  data(myPeakList)
  assignChromosomeRegion(myPeakList, nucleotideLevel=TRUE, precedence=c("Promoters", "Exons", "Introns", "fiveUTRs", "threeUTRs", "immediateDownstream"), TxDb=TxDb)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
