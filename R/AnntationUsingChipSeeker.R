#' AnntationUsingChipSeeker
#' 
#' Use ChipSeeker to annotate data
#'
#' @param dir.name: path for bed files
#' @param input.file.pattern : input file pattern
#' @param out.dir.name output file directory
#' @param txdb Annotation databased to be used("hg19","hg38")
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
#' txdb="hg19"
#' DD=5000
#' AnntationUsingChipSeeker(dir.name,input.file.pattern,out.dir.name,DD)
#' 
#' 
#' dir.name="/Volumes/Bioinformatics$/2016/Danny/Analysis4Peaks/"
#' input.file.pattern="bed"
#' out.dir.name="/Users/axy148/Aimin_project/Danny/"
#' DD=5000
#' 
#' 
#' dir.name="/Volumes/Bioinformatics$/2016/Danny/Analysis4Peaks/common_peaks_bed/"
#' input.file.pattern="bed"
#' out.dir.name="/Volumes/Bioinformatics$/2016/Danny/Analysis4Peaks/"
#' DD=5000
#' AnntationUsingChipSeeker(dir.name,input.file.pattern,out.dir.name,txdb=txdb,DD)
#' 
AnntationUsingChipSeeker <- function(dir.name,input.file.pattern,out.dir.name,txdb=c("hg19","hg38"),DD) {

  re<-ParserReadFiles(dir.name,input.file.pattern,out.dir.name)
  
  re.bed<-re$input
  
  re.peaks.only.bed.2 <- re.bed
  
  # if(length(dir(dir.name,pattern="peaks.bed"))!=0)
  # {
  # re.peaks.only.bed.2<-FL(re.bed,'peaks')
  # cat("peaks\n")
  # print(re.peaks.only.bed.2)
  # }
  # 
  # if(length(dir(dir.name,pattern="summits.bed"))!=0){
  # re.summits.only.bed<-FL(re.bed,'summits')
  # cat("summits\n")
  # print(re.summits.only.bed)
  # }
  
  txdb<-match.arg(txdb)
    
    switch (txdb,
            hg38 = {
              cat("Use hg38\n")
               txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
            },
            {
              cat("Use hg19\n") 
              txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
            }
           )
  
   temp3=file.path(re$output,"Annotation")
   
   if(!dir.exists(temp3)){dir.create(temp3)}
   
    d=DD
   lapply(1:length(re.peaks.only.bed.2),function(u,re.peaks.only.bed.2,d){
     
     peaks=readPeakFile(re.peaks.only.bed.2[[u]],as="data.frame")
     
     print(head(peaks))
     
     peakAnno <- annotatePeak(re.peaks.only.bed.2[[u]], tssRegion=c(-d, d),
                              TxDb=txdb, annoDb="org.Hs.eg.db")
     x_name=names(re.peaks.only.bed.2)[u]
     cat(x_name)
     png(file.path(temp3,paste0(x_name,"_",d,"_around_tss_annotation_pie.png")))             
     plotAnnoPie(peakAnno)
     dev.off()
     
     peaks.anno=as.data.frame(peakAnno)
     
     print(head(peaks.anno))
     
     #print(paste0(peaks[,c(2,3)]))
     print(colnames(peaks.anno))
     write.table(peaks.anno,file=file.path(temp3,paste0(x_name,"_",d,"_around_tss_annotation_4_only_mapped_peaks.xls")),
                 row.names = FALSE,quote=FALSE,sep="\t")
     
     unmapped.peaks<-peaks[-which(paste0(peaks[,2],"_",peaks[,3]) %in% paste0(peaks.anno[,2],"_",peaks.anno[,3])),]
        
     cat(dim(peaks)[1]," ",dim(peaks.anno)[1]," ",dim(unmapped.peaks)[1],"\n")
     
     
     if(dim(unmapped.peaks)[1]!=0){
     
     colnames(unmapped.peaks)=colnames(peaks.anno)[1:6]
     
     unmapped.peaks.3<-smartbind(peaks.anno,unmapped.peaks)
     
     unmapped.peaks.4<-unmapped.peaks.3[order(unmapped.peaks.3[,1],unmapped.peaks.3[,2]),]
     
     write.table(unmapped.peaks.4,file=file.path(temp3,paste0(x_name,"_",d,"_around_tss_annotation_4_all_peaks.xls")),row.names = FALSE,quote=FALSE,sep="\t")
     }
     
     
   },re.peaks.only.bed.2,d)
   
}