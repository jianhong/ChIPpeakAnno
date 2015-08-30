\name{TSS.zebrafish.Zv8}
\Rdversion{1.1}
\alias{TSS.zebrafish.Zv8}
\docType{data}
\title{
  TSS annotation data for zebrafish (Zv8) obtained from biomaRt
}
\description{
 A GRanges object to annotate TSS for zebrafish (Zv8) obtained from biomaRt
}
\usage{data(TSS.zebrafish.Zv8)}
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
                   host="may2009.archive.ensembl.org", 
                   path="/biomart/martservice", 
                   dataset="drerio_gene_ensembl")
   
   getAnnotation(mart, featureType = "TSS")

}
\examples{
data(TSS.zebrafish.Zv8)
slotNames(TSS.zebrafish.Zv8)
}
\keyword{datasets}
