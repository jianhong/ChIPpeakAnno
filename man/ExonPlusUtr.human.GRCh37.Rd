\name{ExonPlusUtr.human.GRCh37}
\alias{ExonPlusUtr.human.GRCh37}
\docType{data}
\title{
 Gene model with exon, 5' UTR and 3' UTR information for human sapiens (GRCh37) 
 obtained from biomaRt
}
\description{
 Gene model with exon, 5' UTR and 3' UTR information for human sapiens (GRCh37) 
 obtained from biomaRt
}
\usage{data(ExonPlusUtr.human.GRCh37)}
\format{
    RangedData with slot start holding the start position of the exon, slot end 
    holding the end position of the exon, slot rownames holding ensembl 
    transcript id and slot space holding the chromosome location where the gene 
    is located. In addition, the following variables are included.
  \describe{
    \item{\code{strand}}{1 for positive strand and -1 for negative strand}
    \item{\code{description}}{description of the transcript}
\item{\code{ensembl_gene_id}}{gene id}
\item{\code{utr5start}}{5' UTR start}
\item{\code{utr5end}}{5' UTR end}
\item{\code{utr3start}}{3' UTR start}
\item{\code{utr3end}}{3' UTR end}}
}
\details{
 used in the examples
  Annotation data obtained by:
  mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  ExonPlusUtr.human.GRCh37 = 
      getAnnotation(mart=human, featureType="ExonPlusUtr")
}
\examples{
data(ExonPlusUtr.human.GRCh37)
slotNames(ExonPlusUtr.human.GRCh37)
}
\keyword{datasets}
