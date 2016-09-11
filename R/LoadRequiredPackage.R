#' Title
#'
#' @return
#' @export
#'
#' @examples
#' LoadRequiredPackage()
#'
LoadRequiredPackage <- function() {
  library(LOLA)

  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  library(clusterProfiler)
  library(regioneR)
  library(BSgenome.Hsapiens.UCSC.hg19)
}
