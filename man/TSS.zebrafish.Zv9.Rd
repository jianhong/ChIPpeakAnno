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
   GRanges with slot start holding the start position of the gene, slot end 
   holding the end position of the gene, slot names holding ensembl gene id,
   slot seqnames holding the chromosome location where the gene is located and 
   slot strand holding the strinad information. In addition, the following 
   variables are included.
  \describe{
    \item{\code{description}}{description of the gene}}
}
\details{
Annotation data obtained by: 

    mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                   host="mar2015.archive.ensembl.org", 
                   path="/biomart/martservice", 
                   dataset="drerio_gene_ensembl")
    
    getAnnotation(mart, featureType = "TSS")
}
\examples{
data(TSS.zebrafish.Zv9)
slotNames(TSS.zebrafish.Zv9)
}
\keyword{datasets}
