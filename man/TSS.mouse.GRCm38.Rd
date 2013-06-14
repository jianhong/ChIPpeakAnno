\name{TSS.mouse.GRCm38}
\Rdversion{1.1}
\alias{TSS.mouse.GRCm38}
\docType{data}
\title{
 TSS annotation data for Mus musculus (GRCm38.p1) obtained from biomaRt
}
\description{
 TSS annotation data for Mus musculus (GRCm38.p1) obtained from biomaRt
}
\usage{data(TSS.mouse.GRCm38)}
\format{
    RangedData with slot start holding the start position of the gene, slot end holding the end position of the gene, slot rownames holding ensembl gene id and slot space holding the chromosome location where the gene is located. In addition, the following variables are included.
  \describe{
     \item{\code{strand}}{1 for positive strand and -1 for negative strand}
    \item{\code{description}}{description of the gene}}}

\details{
 Annotation data obtained by:
    
	 mart = useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

 	getAnnotation(mart, featureType = "TSS")

}
\examples{
data(TSS.mouse.GRCm38)
slotNames(TSS.mouse.GRCm38)
}
\keyword{datasets}
