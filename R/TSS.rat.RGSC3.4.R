#' TSS annotation data for rat (RGSC3.4) obtained from biomaRt
#' 
#' TSS annotation data for rat (RGSC3.4) obtained from biomaRt
#' 
#' Annotation data obtained by:
#' 
#' mart = useMart(biomart = "ensembl", dataset = "rnorvegicus_gene_ensembl")
#' 
#' getAnnotation(mart, featureType = "TSS")
#' 
#' @name TSS.rat.RGSC3.4
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
#' data(TSS.rat.RGSC3.4)
#' slotNames(TSS.rat.RGSC3.4)
"TSS.rat.RGSC3.4"
