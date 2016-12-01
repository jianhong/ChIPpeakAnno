#' AnnotatePeak3
#'
#' @param input.file.dir 
#' @param output.file.dir 
#' @param genome 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' input.file.dir="/media/H_driver/2016/Danny/Danny_chip_PeakCall/"
#' output.file.dir="/media/H_driver/2016/Danny/Danny_chip_PeakCall/"
#'
#' AnnotatePeak3(input.file.dir,output.file.dir,genome="Mm")
#' AnnotatePeak3(input.file.dir,output.file.dir,genome="Hs")
#' 
#' 
#' 
AnnotatePeak3<- function(input.file.dir,output.file.dir,genome) {
  
  library(ChIPpeakAnno)
  
  re<-ParserReadFiles(input.file.dir,"bed",output.file.dir)
  
  #print(file.name.2)
  
  file.name.2<-re$input
  output.dir.name=re$output
  
  #dir.name=input.file.dir
  #input.file.pattern=input.file.pattern
  
  #dir.name=reformatPath(dir.name)
  
  #output.dir.name=reformatPath(output.file.dir)
  
  #print(output.dir.name)
  #temp=Sys.time()
  #temp1=gsub(":","-",Sys.time())
  #temp2=gsub(" ","-",temp1)
  #temp3=paste0(output.dir.name,"AnalysisResults_at_",temp2)
  
  temp3=output.dir.name
  
  dir.create(temp3)
  
  #file.name=paste0(dir.name,dir(dir.name,recursive = TRUE,pattern=input.file.pattern))
  #file.name.2<-as.list(file.name)
  
  #names(file.name.2)=sapply(strsplit(file.name,split="\\/"),"[[",index.file)
  
  #print(file.name.2)
  
  file.name.3<-file.name.2
  
  #sample.name<-sapply(strsplit(names(file.name.3),split="_peaks_"),"[[",1)
  
  #names(file.name.3)=sample.name
  
  file.name.4 <-file.name.3
  
  re.out<-lapply(file.name.4,function(u){
    re=toGRanges(u,format="BED")
    #colnames(re)=c("Count","GeneName")
    re
  })
  
  head(re.out[[1]])
  
  re.out.L<-lapply(re.out,function(u){
    re=length(u)
    #colnames(re)=c("Count","GeneName")
    re
  })
  
  if(genome=="Mm"){
    
    annoData <- toGRanges(EnsDb.Mmusculus.v75, feature="gene")
    
    ol <- findOverlapsOfPeaks(re.out[c(2,4,1)])
    
    overlaps<-ol$peaklist$`11_2470IUPUI_WT_BM_SMC1_peaks.bed///13_2470IUPUI_WT_BM_Rad21_peaks.bed///10_WT_BM_ASXL1_peaks.bed`
    
    binOverFeature(overlaps, annotationData=annoData,
                   radius=5000, nbins=20, FUN=length, errFun=0,
                   ylab="count",
                   main="Distribution of aggregated peak numbers around TSS")
    
    overlaps.trimmed<-trim(overlaps, use.names=TRUE)
    
    library(EnsDb.Mmusculus.v79)
    dd.GRCm39.mm10<-toGRanges(EnsDb.Mmusculus.v75)
    #seqinfo(dd.GRCm39.mm10)
    #seqlevels(dd.GRCm39.mm10)
    
    seqlevels(dd.GRCm39.mm10,force=TRUE) <- c("chr1","chr10","chr11","chr12","chr13",
                                              "chr14","chr15","chr16","chr17","chr18","chr19","chr2",
                                              "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY")
    
    seqinfo(overlaps)<-seqinfo(dd.GRCm39.mm10)
    #GRCm38/mm10
    #dd<-toGRanges(EnsDb.Mmusculus.v79)
    #seqinfo(dd)
    #library(ensembldb)
    #library(GenomeInfoDb)
    seqlevelsStyle(overlaps.trimmed) <- seqlevelsStyle(dd.GRCm39.mm10)
    overlaps.anno<-annoPeaks(overlaps.trimmed,dd.GRCm39.mm10)
    write.table(overlaps.anno,file=paste0(temp3,"/","annotation.txt"),row.names = FALSE,quote=FALSE,sep="\t")
  }else if(genome=="Hs"){
    
    library(EnsDb.Hsapiens.v75)
    #annoData<-toGRanges(EnsDb.Hsapiens.v75, feature="gene")
    
    dd.hs<-toGRanges(EnsDb.Hsapiens.v75)
    
    print(seqinfo(dd.hs))
    print(seqlevels(dd.hs))
    
    #print(seqlevels(dd.hs)[,1])
    
    #print(seqlevels(re.out[[1]])[,1])
    
    # seqlevels(dd.hs,force=TRUE) <- c("chr1","chr10","chr11","chr12","chr13",
    #                                           "chr14","chr15","chr16","chr17","chr18","chr19","chr2",
    #                                           "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY")
    
    #temp4=
    
    re.out.L<-lapply(1:length(re.out),function(u,re.out,dd.hs){
      
      x=re.out[[u]]
      x_name=names(re.out)[u]
      
      seqlevels(dd.hs,force=TRUE)<-seqinfo(x)@seqnames
      #print(seqinfo(re.out.trimmed))
      #print(seqlevels(re.out.trimmed))
      seqinfo(x)<-seqinfo(dd.hs)
      #GRCm38/mm10
      #dd<-toGRanges(EnsDb.Mmusculus.v79)
      #seqinfo(dd)
      #library(ensembldb)
      #library(GenomeInfoDb)
      seqlevelsStyle(x) <- seqlevelsStyle(dd.hs)
      re.out.trimmed<-trim(x, use.names=TRUE)
      overlaps.anno<-annoPeaks(re.out.trimmed,dd.hs)
      
      write.table(overlaps.anno,file=paste0(temp3,"/",x_name,"_annotation.txt"),row.names = FALSE,quote=FALSE,sep="\t")
    },re.out,dd.hs)
    
  }
  
}