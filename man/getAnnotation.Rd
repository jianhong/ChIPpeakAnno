\name{getAnnotation}
\alias{getAnnotation}
\title{ Obtain the TSS, exon or miRNA annotation for the specified species }
\description{
  Obtain the TSS, exon or miRNA annotation for the specified species using the biomaRt package
}
\usage{
getAnnotation(mart, 
featureType=c("TSS","miRNA", "Exon", "5utr", "3utr", "ExonPlusUtr", "transcript"), 
output=c("GRanges", "RangedData"))
}
\arguments{
  \item{mart}{A mart object, see useMart of biomaRt package for details. }
  \item{featureType}{TSS, miRNA, Exon, 5'UTR, 3'UTR, transcript or Exon plus UTR. The default is TSS. }
  \item{output}{The class of the output data, \link[GenomicRanges:GRanges]{GRanges} or \link[IRanges:RangedData]{RangedData}. The default is RangedData.}
}
\value{
  \link[GenomicRanges:GRanges]{GRanges} or \link[IRanges:RangedData]{RangedData} with slot start holding the start position of the feature, slot end holding the end position of the feature, slot names holding the id of the feature, slot space holding the chromosome location where the feature is located. In addition, the following variables are included.
    \item{\code{strand}}{1 for positive strand and -1 for negative strand where the feature is located}
    \item{\code{description}}{description of the feeature such as gene}
}

\references{ Durinck S. et al. (2005) BioMart and Bioconductor: a powerful link between biological biomarts and microarray data analysis. Bioinformatics, 21, 3439-3440.  }

\author{ Lihua Julie Zhu, Jianhong Ou }
\note{ 
	For featureType of TSS, start is the transcription start site if strand is 1 (plus strand), otherwise, end is the transcription start site
}
\examples{
if (interactive())
{
	mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
	Annotation <- getAnnotation(mart, featureType="TSS", output="GRanges")
}
}
\keyword{ misc }
