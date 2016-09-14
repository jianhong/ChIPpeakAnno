#' InstallRequiredPackage
#'
#' @return
#' @export
#'
#' @examples
#' InstallRequiredPackage()
#'
InstallRequiredPackage <- function() {
  source("http://bioconductor.org/biocLite.R")
  biocLite("GenomicRanges")
  devtools::install_github("nsheff/LOLA")
  biocLite("ChIPseeker")
  biocLite("regioneR")
  biocLite("DiffBind")
  library(BiocInstaller)
  biocLite("BSgenome.Hsapiens.UCSC.hg19")
  library(DiffBind)
}
