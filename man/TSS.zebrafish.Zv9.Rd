\name{TSS.zebrafish.Zv9}
\alias{TSS.zebrafish.Zv9}
\docType{data}
\title{
TSS annotation for Danio rerio (Zv9) obtained from biomaRt
}
\description{
TSS annotation for Danio rerio (Zv9) obtained from biomaRt
}
\usage{data(TSS.zebrafish.Zv9)}
\format{
   RangedData with slot start holding the start position of the gene, slot end holding the end position of the gene, slot rownames holding ensembl gene id and slot space holding the chromosome location where the gene is located. In addition, the following variables are included.
  \describe{
    \item{\code{strand}}{1 for positive strand and -1 for negative strand}
    \item{\code{description}}{description of the gene}}
}
\details{
Annotation data obtained by: 

    mart = useMart(biomart = "ensembl", dataset = "drerio_gene_ensembl")
	
	getAnnotation(mart, featureType = "TSS")
}
\examples{
data(TSS.zebrafish.Zv9)
slotNames(TSS.zebrafish.Zv9)
}
\keyword{datasets}
