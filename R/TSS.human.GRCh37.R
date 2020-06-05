#' TSS annotation for human sapiens (GRCh37) obtained from biomaRt
#' 
#' TSS annotation for human sapiens (GRCh37) obtained from biomaRt
#' 
#' The dataset TSS.human.GRCh37 was obtained by:
#' 
#' mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",
#' path="/biomart/martservice", dataset = "hsapiens_gene_ensembl")
#' 
#' getAnnotation(mart, featureType = "TSS")
#' 
#' @name TSS.human.GRCh37
#' @docType data
#' @format A GRanges object with slot start holding the start position of the
#' gene, slot end holding the end position of the gene, slot names holding
#' ensembl gene id, slot seqnames holding the chromosome location where the
#' gene is located and slot strand holding the strinad information. In
#' addition, the following variables are included.  \describe{
#' \item{list("description")}{description of the gene}}
#' @keywords datasets
#' @examples
#' 
#' data(TSS.human.GRCh37)
#' slotNames(TSS.human.GRCh37)
#' 
"TSS.human.GRCh37"