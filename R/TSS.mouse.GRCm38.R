#' TSS annotation data for Mus musculus (GRCm38.p1) obtained from biomaRt
#' 
#' TSS annotation data for Mus musculus (GRCm38.p1) obtained from biomaRt
#' 
#' Annotation data obtained by:
#' 
#' mart = useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
#' 
#' getAnnotation(mart, featureType = "TSS")
#' 
#' @name TSS.mouse.GRCm38
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
#' data(TSS.mouse.GRCm38)
#' slotNames(TSS.mouse.GRCm38)
#' 
"TSS.mouse.GRCm38"
