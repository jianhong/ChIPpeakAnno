#' TSS annotation for human sapiens (GRCh38) obtained from biomaRt
#' 
#' TSS annotation for human sapiens (GRCh38) obtained from biomaRt
#' 
#' used in the examples Annotation data obtained by:
#' 
#' mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#' 
#' getAnnotation(mart, featureType = "TSS")
#' 
#' @name TSS.human.GRCh38
#' @docType data
#' @format A 'GRanges' [package "GenomicRanges"] object with ensembl id as
#' names.
#' @keywords datasets
#' @examples
#' 
#' data(TSS.human.GRCh38)
#' slotNames(TSS.human.GRCh38)
#' 
"TSS.human.GRCh38"
