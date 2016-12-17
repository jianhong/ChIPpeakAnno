#'  GetResultsFromDiffBind
#'  
#'  Use DiffBind to process peak profiles
#'
#' @param mcf7
#' @param output.file.dir 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' data(mcf7)
#' output.file.dir="/Volumes/Bioinformatics$/2016/Danny/Analysis4Peaks/"
#' GetResultsFromDiffBind(mcf7,output.file.dir)
#' 
#' 
GetResultsFromDiffBind<-function(mcf7,output.file.dir){
  
  temp<-mcf7
  
  sampID.v<-colnames(temp$class) 
  sampID.v.2<-unlist(lapply(1:length(sampID.v),function(u,sampID.v){
    
    x=sampID.v[u]
    
    pos=regexpr("_",x)
    
    y=substr(x,27,pos-1)
    y
  },sampID.v))
  
  #sampID.v.2<-substr(sampID.v,27,nchar(sampID.v)-3)
 
  colnames(temp$class)<-sampID.v.2
  colnames(temp$class)<-sampID.v.2
  temp$class[1,]<-sampID.v.2
  
  #Check the correlation between peaks
  png(paste0(output.file.dir,"CorrelationHeatmap.png"))
  dba.plotHeatmap(temp,margin=15)
  dev.off()
  
  #Merge the replicates of each set
  temp2<-dba.peakset(temp,consensus = -DBA_REPLICATE)
  temp22<-dba(temp2,mask =c(9,16,17:23))
  
  Tissue1<-temp22$class[2,]
  Tissue2<-unique(temp22$class[2,])
  
  TF<-unique(temp22$class[3,])
  TF.n<-length(TF)
  
  #Draw 
  png(paste0(output.file.dir,"Venn.png"))
  
  par(mfrow=c(2,2))
  for(i in 1:length(Tissue2)){
    po<-which(Tissue1 %in% Tissue2[i])
    if(length(po)>1)
    {
      dba.plotVenn(temp22,mask=po,main=Tissue2[i])
    }else
    {
      cat(paste0("For ",Tissue2[i],": Only one peak profile for ",TF.n," TFs"))
    }
  }
  dev.off()
  #Get common peaks for TFs
  
  p.common<-lapply(1:length(Tissue2),function(u,Tissue1,Tissue2,temp22){
    
    po<-which(Tissue1 %in% Tissue2[u])
    
    if(length(po)>1)
    {
      common.peaks<-dba.overlap(temp22,mask=po)
      y<-common.peaks$inAll
    }else{
      y<-NULL}
    y
  },Tissue1,Tissue2,temp22)
  
  names(p.common)<-Tissue2
  
  p.common<-p.common[lapply(p.common,length) > 0]
  
  library(EnsDb.Hsapiens.v75)
  dd.hs<-toGRanges(EnsDb.Hsapiens.v75)
  
  lapply(1:length(p.common),function(u,p.common,dd.hs,output.file.dir){
    x=p.common[[u]]
    x_name=names(p.common)[u]
    seqlevels(dd.hs,force=TRUE)<-seqinfo(x)@seqnames
    seqlevelsStyle(x) <- seqlevelsStyle(dd.hs)
    re.out.trimmed<-trim(x, use.names=TRUE)
    overlaps.anno<-annoPeaks(re.out.trimmed,dd.hs)
    
    write.table(overlaps.anno,file=paste0(output.file.dir,x_name,"_annotation_2.txt"),row.names = FALSE,quote=FALSE,sep="\t")
  },p.common,dd.hs,output.file.dir)
  

  # temp3<-dba.peakset(temp2,bRetrieve = TRUE)
  # temp.no.replicate<-dba(dba.peakset(temp2,c(9,16),bRetrieve=TRUE,DataTyp=DBA_DATA_FRAME))
  # temp.replicate<-dba.peakset(temp2,17:23,bRetrieve=TRUE,DataTyp=DBA_DATA_FRAME)
  # temp.c<-dba.peakset(peaks = temp.replicate)
  # temp3<-dba.peakset(temp2,temp2$masks$Consensus,bRetrieve = TRUE)
  # temp3<-dba.peakset(NULL,peaks=temp2$masks$Consensus)
  # 
}