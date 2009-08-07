\name{TSS.human.NCBI36}
\Rdversion{1.1}
\alias{TSS.human.NCBI36}
\docType{data}
\title{
  TSS annotation for human sapiens (NCBI36) obtained from biomaRt
}
\description{
  TSS annotation for human sapiens (NCBI36) obtained from biomaRt
}
\usage{data(TSS.human.NCBI36)}
\format{
   RangedData with slot start holding the start position of the gene, slot end holding the end position of the gene, slot rownames holding ensembl gene id and slot space holding the chromosome location where the gene is located. In addition, the following variables are included.
  \describe{
    \item{\code{strand}}{1 for positive strand and -1 for negative strand}
    \item{\code{description}}{description of the gene}}}
\details{
 used in the examples
  Annotation data obtained by:

	mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
		
	getAnnotation(mart, featureType = "TSS")

}
\source{

}
\references{
}
\examples{
data(TSS.human.NCBI36)
slotNames(TSS.human.NCBI36)
}
\keyword{datasets}
