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
}
