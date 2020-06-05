#' TSS annotation for Danio rerio (Zv9) obtained from biomaRt
#' 
#' TSS annotation for Danio rerio (Zv9) obtained from biomaRt
#' 
#' Annotation data obtained by:
#' 
#' mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
#' host="mar2015.archive.ensembl.org", path="/biomart/martservice",
#' dataset="drerio_gene_ensembl")
#' 
#' getAnnotation(mart, featureType = "TSS")
#' 
#' @name TSS.zebrafish.Zv9
#' @docType data
#' @format GRanges with slot start holding the start position of the gene, slot
#' end holding the end position of the gene, slot names holding ensembl gene
#' id, slot seqnames holding the chromosome location where the gene is located
#' and slot strand holding the strinad information. In addition, the following
#' variables are included.  \describe{ \item{list("description")}{description
#' of the gene}}
#' @keywords datasets
#' @examples
#' 
#' data(TSS.zebrafish.Zv9)
#' slotNames(TSS.zebrafish.Zv9)
#' 
"TSS.zebrafish.Zv9"
