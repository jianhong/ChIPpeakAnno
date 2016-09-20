#' Title
#'
#' @return
#' @export
#'
#' @examples
#'
#' which.bed.file="/media/H_driver/2016/Yang/Results/WT_triple_overlap.bed"
#'
ProcessUsingChipSeeker <- function(which.bed.file) {

  #Example

  files <- getSampleFiles()
  print(files)

  peak <- readPeakFile(files[[4]])
  peak


  data("tagMatrixList")
  tagMatrix <- tagMatrixList[[4]]

  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

  tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")

  peakAnno <- annotatePeak(files[[4]], tssRegion=c(-3000, 3000),
                           TxDb=txdb, annoDb="org.Hs.eg.db")

  plotAnnoPie(peakAnno)
  peak

  #Yang data
  file.triple=which.bed.file

  peak.3<-readPeakFile(file.triple)
  txdb.mm9 <- TxDb.Mmusculus.UCSC.mm9.knownGene

  peakAnno.triple <- annotatePeak(file.triple, tssRegion=c(-1000, 1000),
                           TxDb=txdb.mm9)

  plotAnnoPie(peakAnno.triple)

  getGEOspecies()
  getGEOgenomeVersion()

  mm9 <- getGEOInfo(genome="mm9", simplify=TRUE)

  head(mm9)
  gsm <- mm9$gsm[sample(nrow(mm9), 10)]

  downloadGEObedFiles(genome=gsm, destDir="mm9")

}
