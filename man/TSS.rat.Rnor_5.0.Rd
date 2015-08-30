\name{TSS.rat.Rnor_5.0}
\Rdversion{1.1}
\alias{TSS.rat.Rnor_5.0}
\docType{data}
\title{
 TSS annotation data for Rattus norvegicus (Rnor_5.0) obtained from biomaRt
}
\description{
 TSS annotation data for Rattus norvegicus (Rnor_5.0) obtained from biomaRt
}
\usage{data(TSS.rat.Rnor_5.0)}
\format{
  GRanges with slot start holding the start position of the gene, slot end 
  holding the end position of the gene, slot names holding ensembl gene id, 
  slot seqnames holding the chromosome location where the gene is located and 
  slot strand holding the strinad information. In addition, the following 
  variables are included.
  \describe{
    \item{\code{description}}{description of the gene}}}

\details{
 Annotation data obtained by:
 
    mart = useMart(biomart = "ensembl", dataset = "rnorvegicus_gene_ensembl")

    getAnnotation(mart, featureType = "TSS")
}
\examples{
data(TSS.rat.Rnor_5.0)
slotNames(TSS.rat.Rnor_5.0)
}
\keyword{datasets}
