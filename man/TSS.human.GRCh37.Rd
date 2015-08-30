\name{TSS.human.GRCh37}
\alias{TSS.human.GRCh37}
\docType{data}
\title{
  TSS annotation for human sapiens (GRCh37) obtained from biomaRt
}
\description{
  TSS annotation for human sapiens (GRCh37) obtained from biomaRt
}
\usage{data(TSS.human.GRCh37)}
\format{
   A GRanges object with slot start holding the start position of the gene, 
   slot end holding the end position of the gene, slot names holding ensembl 
   gene id, slot seqnames holding the chromosome location where the gene is 
   located and slot strand holding the strinad information. In addition, the 
   following variables are included.
  \describe{
    \item{\code{description}}{description of the gene}}}
\details{
 The dataset TSS.human.GRCh37 was obtained by:

    mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",
        path="/biomart/martservice", dataset = "hsapiens_gene_ensembl")
        
    getAnnotation(mart, featureType = "TSS")

}
\examples{
data(TSS.human.GRCh37)
slotNames(TSS.human.GRCh37)
}
\keyword{datasets}
