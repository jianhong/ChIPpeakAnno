#' TSS annotation for human sapiens (NCBI36) obtained from biomaRt
#' 
#' TSS annotation for human sapiens (NCBI36) obtained from biomaRt
#' 
#' used in the examples Annotation data obtained by:
#' 
#' mart = useMart(biomart = "ensembl_mart_47", dataset =
#' "hsapiens_gene_ensembl", archive=TRUE)
#' 
#' getAnnotation(mart, featureType = "TSS")
#' 
#' @name TSS.human.NCBI36
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
#' data(TSS.human.NCBI36)
#' slotNames(TSS.human.NCBI36)
#' 
"TSS.human.NCBI36"
