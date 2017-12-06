\name{assignChromosomeRegion}
\alias{assignChromosomeRegion}
\title{
Summarize peak distribution over exon, intron, enhancer, proximal promoter, 
5 prime UTR and 3 prime UTR
}
\description{
Summarize peak distribution over exon, intron, enhancer, proximal promoter, 
5 prime UTR and 3 prime UTR
}
\usage{
assignChromosomeRegion(peaks.RD, exon, TSS, utr5, utr3, 
         proximal.promoter.cutoff=1000L, immediate.downstream.cutoff=1000L, 
         nucleotideLevel=FALSE, precedence=NULL, TxDb=NULL)
}
\arguments{
  \item{peaks.RD}{
peaks in GRanges: See example below
}
  \item{exon}{
exon data obtained from getAnnotation or customized annotation of class GRanges 
containing additional variable:  strand (1 or + for plus strand and -1 or - for
minus strand). This parameter is for backward compatibility only. 
\code{\link[GenomicFeatures]{TxDb}} should be used instead.
}
  \item{TSS}{
TSS data obtained from getAnnotation or customized annotation of class GRanges 
containing additional variable:  strand (1 or + for plus strand and -1 or - for
minus strand). For example, data(TSS.human.NCBI36),data(TSS.mouse.NCBIM37), 
data(TSS.rat.RGSC3.4) and data(TSS.zebrafish.Zv8). This parameter is for 
backward compatibility only. \code{\link[GenomicFeatures]{TxDb}} should be 
used instead.
}
  \item{utr5}{
5 prime UTR data obtained from getAnnotation or customized annotation of class 
GRanges containing additional variable:  strand (1 or + for plus strand and -1
or - for minus strand). This parameter is for backward compatibility only. 
\code{\link[GenomicFeatures]{TxDb}} should be used instead.
}
  \item{utr3}{
3 prime UTR data obtained from getAnnotation or customized annotation of class
GRanges containing additional variable:  strand (1 or + for plus strand and -1 
or - for minus strand). This parameter is for backward compatibility only. 
\code{\link[GenomicFeatures]{TxDb}} should be used instead.
}
  \item{proximal.promoter.cutoff}{
Specify the cutoff in bases to classify proximal promoter or enhencer. Peaks 
that reside within proximal.promoter.cutoff upstream from or overlap with 
transcription start site are classified as proximal promoters. Peaks that 
reside upstream of the proximal.promoter.cutoff from gene start are classified
as enhancers. The default is 1000 bases.

}
  \item{immediate.downstream.cutoff}{
Specify the cutoff in bases to classify immediate downstream region or enhancer
region. Peaks that reside within immediate.downstream.cutoff downstream of 
gene end but not overlap 3 prime UTR are classified as immediate downstream. 
Peaks that reside downstream over immediate.downstreatm.cutoff from gene end 
are classified as enhancers. The default is 1000 bases.
}
  \item{nucleotideLevel}{Logical. Choose between peak centric and nucleotide 
  centric view. Default=FALSE}
  \item{precedence}{If no precedence specified,
  double count will be enabled, which means that if a peak overlap with both
promoter and 5'UTR, both promoter and 5'UTR will be incremented. If a
precedence order is specified, for example, if promoter is specified before
5'UTR, then only promoter will be incremented for the same example. 
The values could be any conbinations of "Promoters", "immediateDownstream", 
"fiveUTRs", "threeUTRs", "Exons" and "Introns", Default=NULL}
  \item{TxDb}{an object of \code{\link[GenomicFeatures]{TxDb}}}
}
\value{
A list of two named vectors: percentage and jaccard (Jaccard Index). The 
information in the vectors:  
 \item{\code{Exons}}{Percent or the picard index of the peaks resided in 
 exon regions.}
    \item{\code{Introns}}{Percent or the picard index of the peaks resided in 
    intron regions.}
    \item{\code{fiveUTRs}}{Percent or the picard index of the peaks resided in 
    5 prime UTR regions.}
    \item{\code{threeUTRs}}{Percent or the picard index of the peaks resided 
    in 3 prime UTR regions.}
   \item{\code{Promoter}}{Percent or the picard index of the peaks resided in 
   proximal promoter regions.}
   \item{\code{ImmediateDownstream}}{Percent or the picard index of the peaks 
   resided in immediate downstream regions.}
   \item{\code{Enhancer.Silencer}}{Percent or the picard index of the peaks 
   resided in enhancer/silencer regions.}
}
\references{ 1. Zhu L.J. et al. (2010) ChIPpeakAnno: a Bioconductor package to
annotate ChIP-seq and ChIP-chip data. BMC Bioinformatics 2010, 
11:237doi:10.1186/1471-2105-11-237

2. Zhu L.J. (2013) Integrative analysis of ChIP-chip and ChIP-seq dataset. 
Methods Mol Biol. 2013;1067:105-24. doi: 10.1007/978-1-62703-607-8\_8.}
\author{ Jianhong Ou, Lihua Julie Zhu }

\seealso{annotatePeakInBatch, findOverlapsOfPeaks,getEnriched, 
makeVennDiagram,addGeneIDs, peaksNearBDP,summarizePatternInPeaks}

\examples{
if (Sys.getenv("USER")=="jianhongou"){
    ##Display the list of genomes available at UCSC:
    #library(rtracklayer)
    #ucscGenomes()[, "db"]
    ## Display the list of Tracks supported by makeTxDbFromUCSC()
    #supportedUCSCtables()
    ##Retrieving a full transcript dataset for Human from UCSC
    ##TranscriptDb <- 
    ##     makeTxDbFromUCSC(genome="hg19", tablename="ensGene")
    if(require(TxDb.Hsapiens.UCSC.hg19.knownGene)){
        TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene
        exons <- exons(TxDb, columns=NULL)
        fiveUTRs <- unique(unlist(fiveUTRsByTranscript(TxDb)))
        Feature.distribution <- 
            assignChromosomeRegion(exons, nucleotideLevel=TRUE, TxDb=TxDb)
        barplot(Feature.distribution$percentage)
        assignChromosomeRegion(fiveUTRs, nucleotideLevel=FALSE, TxDb=TxDb)
        data(myPeakList)
        assignChromosomeRegion(myPeakList, nucleotideLevel=TRUE, 
                               precedence=c("Promoters", "immediateDownstream", 
                                            "fiveUTRs", "threeUTRs", 
                                            "Exons", "Introns"), 
                               TxDb=TxDb)
    }
}
}
\keyword{ misc }
