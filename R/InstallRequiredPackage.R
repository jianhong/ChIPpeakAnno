#' InstallRequiredPackage
#'
#' @return
#' @export
#'
#' @examples
#' InstallRequiredPackage()
#'
InstallRequiredPackage <- function() {

  #install from bioc
  
  source("http://bioconductor.org/biocLite.R")
  biocLite("GenomicRanges")
  biocLite("ChIPseeker")
  biocLite("regioneR")
  biocLite("DiffBind")
  biocLite("BSgenome.Hsapiens.UCSC.hg19")
  biocLite("TxDb.Mmusculus.UCSC.mm9.knownGene")
  biocLite("EnsDb.Hsapiens.v75")
  devtools::install_github("sheffien/simpleCache")
  biocLite("LOLA")
  biocLite("ChIPpeakAnno")
  biocLite("BSgenome.Mmusculus.UCSC.mm9")
  biocLite("motifStack")
  biocLite("BSgenome.Ecoli.NCBI.20080805")
  biocLite("EnsDb.Mmusculus.v75")
  biocLite("EnsDb.Mmusculus.v79")
  biocLite("BSgenome.Mmusculus.UCSC.mm10")
  biocLite("org.Mm.eg.db")
  
  #install from git
  devtools::install_github("nsheff/LOLA")
  
  
  #loading
  library(BiocInstaller)
  library(BiocInstaller)
  library(ChIPseeker)
  library(DiffBind)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(TxDb.Mmusculus.UCSC.mm9.knownGene)
  library(LOLA)
  library("simpleCache")
  library(ChIPpeakAnno)
  library(rtracklayer)
  library(org.Mm.eg.db)
  library(motifStack)

  library(BSgenome.Mmusculus.UCSC.mm10)
  library(BSgenome.Mmusculus.UCSC.mm9)
  library(EnsDb.Mmusculus.v75)
  library(EnsDb.Mmusculus.v79)
  
  library(ggplot2)
}
