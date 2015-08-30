\name{getAllPeakSequence}
\alias{getAllPeakSequence}
\title{ Obtain genomic sequences around the peaks
 }
\description{
  Obtain genomic sequences around the peaks leveraging the BSgenome and biomaRt package
}
\usage{
getAllPeakSequence(myPeakList, upstream = 200L, downstream = upstream, 
                   genome, AnnotationData)
}
\arguments{
  \item{myPeakList}{ An object of \link[GenomicRanges:GRanges]{GRanges}: 
  See example below}
  \item{upstream}{upstream offset from the peak start, e.g., 200 }
  \item{downstream}{ downstream offset from the peak end, e.g., 200 }
  \item{genome}{ BSgenome object or mart object. Please refer to 
  available.genomes in BSgenome package and useMart in bioMaRt package 
  for details}
   \item{AnnotationData}{ RangedData used if mart object is parsed in which
   can be obtained from getAnnotation with featureType="TSS". For example, 
   data(TSS.human.NCBI36), data(TSS.mouse.NCBIM37), data(GO.rat.RGSC3.4) and 
   data(TSS.zebrafish.Zv8).  If not supplied, then annotation will be
   obtained from biomaRt automatically using the mart object}
}
\value{
 \link[GenomicRanges:GRanges]{GRanges} with slot start holding the start 
 position of the peak, slot end holding the end position of the peak, 
 slot rownames holding the id of the peak and slot seqnames holding the 
 chromosome where the peak is located. In addition, the following variables 
 are included:
 \item{upstream}{upstream offset from the peak start}
 \item{downstream}{downstream offset from the peak end}
\item{sequence}{the sequence obtained}
}
\references{ Durinck S. et al. (2005) BioMart and Bioconductor: a powerful 
link between biological biomarts and microarray data analysis. 
Bioinformatics, 21, 3439-3440. }
\author{ Lihua Julie Zhu, Jianhong Ou}
\examples{
#### use Annotation data from BSgenome
peaks <- GRanges(seqnames=c("NC_008253", "NC_010468"),
                 IRanges(start=c(100, 500), end=c(300, 600), 
                         names=c("peak1", "peak2")))
library(BSgenome.Ecoli.NCBI.20080805)
seq <- getAllPeakSequence(peaks, upstream=20, downstream=20, genome=Ecoli)
write2FASTA(seq, file="test.fa")
}
\keyword{ misc }
