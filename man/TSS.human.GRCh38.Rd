\name{TSS.human.GRCh38}
\alias{TSS.human.GRCh38}
\docType{data}
\title{
  TSS annotation for human sapiens (GRCh38) obtained from biomaRt
}
\description{
  TSS annotation for human sapiens (GRCh38) obtained from biomaRt
}
\usage{data(TSS.human.GRCh38)}
\format{
  A 'GRanges' [package "GenomicRanges"] object with ensembl id as names.
}
\details{
 used in the examples
  Annotation data obtained by:

    mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
        
    getAnnotation(mart, featureType = "TSS")

}
\examples{
data(TSS.human.GRCh38)
slotNames(TSS.human.GRCh38)
}
\keyword{datasets}
