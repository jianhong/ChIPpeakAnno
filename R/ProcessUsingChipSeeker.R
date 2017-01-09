#' AnntationUsingChipSeeker
#' 
#' Use ChipSeeker to annotate data
#'
#' @param dir.name: path for bed files
#' @param input.file.pattern : input file pattern
#' @param out.dir.name output file directory
#' @param DD distance definition around TSS
#' 
#' @return
#' @export
#'
#' @examples
#'
#' which.bed.file="/Users/axy148/Aimin_project/Danny/bed"
#' 
#' dir.name="/Users/axy148/Aimin_project/Danny/bed/"
#' input.file.pattern="bed"
#' out.dir.name="/Users/axy148/Aimin_project/Danny/"
#' DD=5000
#' AnntationUsingChipSeeker(dir.name,input.file.pattern,out.dir.name,DD)
#' 
#' 
#' dir.name="/Volumes/Bioinformatics$/2016/Danny/Analysis4Peaks/"
#' input.file.pattern="bed"
#' out.dir.name="/Users/axy148/Aimin_project/Danny/"
#' DD=5000
#' 
#' AnntationUsingChipSeeker(dir.name,input.file.pattern,out.dir.name,DD)
#' 
#'
AnntationUsingChipSeeker <- function(dir.name,input.file.pattern,out.dir.name,DD) {

  re<-ParserReadFiles(dir.name,input.file.pattern,out.dir.name)
  
  re.bed<-re$input
  
  if(length(dir(dir.name,pattern="peaks.bed"))!=0)
  {
  re.peaks.only.bed.2<-FL(re.bed,'peaks')
  cat("peaks\n")
  print(re.peaks.only.bed.2)
  }
  
  
  if(length(dir(dir.name,pattern="summits.bed"))!=0){
  re.summits.only.bed<-FL(re.bed,'summits')
  cat("summits\n")
  print(re.summits.only.bed)
  }
  
  # 
  # 
  # 
  # #Example
  # 
  # files <- getSampleFiles()
  # print(files)
  
   peak <- readPeakFile(re.peaks.only.bed.2[[1]])
  # peak
  # 
  # 
  # data("tagMatrixList")
  # tagMatrix <- tagMatrixList[[4]]
  # 
   txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  # 
  # tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
  # 
   
   temp3=paste0(re$output,"ChIPseeker")
   dir.create(temp3)
   
    d=DD
   lapply(1:length(re.peaks.only.bed.2),function(u,re.peaks.only.bed.2,d){
     
     peaks=readPeakFile(re.peaks.only.bed.2[[u]],as="data.frame")
     
     print(head(peaks))
     
     peakAnno <- annotatePeak(re.peaks.only.bed.2[[u]], tssRegion=c(-d, d),
                              TxDb=txdb, annoDb="org.Hs.eg.db")
     x_name=names(re.peaks.only.bed.2)[u]
     cat(x_name)
     png(paste0(temp3,"/",x_name,"_",d,"_around_tss_annotation_pie.png"))             
     plotAnnoPie(peakAnno)
     dev.off()
     
     peaks.anno=as.data.frame(peakAnno)
     
     print(head(peaks.anno))
     
     #print(paste0(peaks[,c(2,3)]))
     
     unmapped.peaks<-peaks[-which(paste0(peaks[,2],"_",peaks[,3]) %in% paste0(peaks.anno[,2],"_",peaks.anno[,3])),]
        
     cat(dim(peaks)[1]," ",dim(peaks.anno)[1]," ",dim(unmapped.peaks)[1],"\n")
     
     
     #unmapped.peaks.part<-as.data.frame(matrix(rep("NA",dim(unmapped.peaks)[1]*14),dim(unmapped.peaks)[1],14))
     
     #unmapped.peaks.2<-cbind(unmapped.peaks,unmapped.peaks.part)
     
     colnames(unmapped.peaks)=colnames(peaks.anno)[1:6]
     
     #unmapped.peaks.3<-rbind.data.frame(unmapped.peaks.2,peaks.anno,stringsAsFactors = FALSE)
     
     unmapped.peaks.3<-smartbind(peaks.anno,unmapped.peaks)
     
     unmapped.peaks.4<-unmapped.peaks.3[order(unmapped.peaks.3[,1],unmapped.peaks.3[,2]),]
     
     write.table(peaks.anno,file=paste0(temp3,"/",x_name,"_",d,"_around_tss_annotation_4_only_mapped_peaks.xls"),row.names = FALSE,quote=FALSE,sep="\t")
     
     write.table(unmapped.peaks.4,file=paste0(temp3,"/",x_name,"_",d,"_around_tss_annotation_4_all_peaks.xls"),row.names = FALSE,quote=FALSE,sep="\t")
     
   },re.peaks.only.bed.2,d)
   
   
   
   #to convert csAnno to GRanges instance, and as.data.frame to convert csAnno to data.frame which can be exported to file by write.table.
   
  # peak
  # 
  # #Yang data
  # file.triple=which.bed.file
  # 
  # peak.3<-readPeakFile(file.triple)
  # txdb.mm9 <- TxDb.Mmusculus.UCSC.mm9.knownGene
  # 
  # peakAnno.triple <- annotatePeak(file.triple, tssRegion=c(-1000, 1000),
  #                          TxDb=txdb.mm9)
  # 
  # plotAnnoPie(peakAnno.triple)
  # 
  # getGEOspecies()
  # getGEOgenomeVersion()
  # 
  # mm9 <- getGEOInfo(genome="mm9", simplify=TRUE)
  # 
  # head(mm9)
  # gsm <- mm9$gsm[sample(nrow(mm9), 10)]
  # 
  # downloadGEObedFiles(genome=gsm, destDir="mm9")

}
