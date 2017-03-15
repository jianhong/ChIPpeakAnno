#'  GetResultsFromDiffBind
#'  
#'  Use DiffBind to process peak profiles
#'
#' @param mcf7
#' @param Mergereplicates
#' @param output.file.dir 
#'
#' @return
#' @export
#'
#' @examples
#' output.file.dir="/Volumes/Bioinformatics$/2017/DannyNewData/BindDiff"
#' mcf7=resmcf
#' re<-GetResultsFromDiffBind(mcf7,"yes",output.file.dir)
#' 
GetResultsFromDiffBind<-function(mcf7,Mergereplicates=c("yes","no"),output.file.dir){
  
  temp<-mcf7
  
  sampID.v<-colnames(temp$class) 
  sampID.v.2<-unlist(lapply(1:length(sampID.v),function(u,sampID.v){
    
    x=sampID.v[u]
    
    y=x
    y
  },sampID.v))
  

  colnames(temp$class)<-sampID.v.2
  temp$class[1,]<-sampID.v.2
  
  temp<-dba(temp,mask =c(1,2,13,14,11,12,15,16))
  
  #Merge the replicates of each set 
  if(Mergereplicates=="yes"){
  temp2<-dba.peakset(temp,consensus = -DBA_REPLICATE)
  #temp22<-dba(temp2,mask =c(9,16,17:23))
  print(temp2)
  ##need to set the flexiable number identified from data sets
  temp22<-dba(temp2,mask =c(9:12))
  }else{
    temp22<-temp
  }
  
  print(temp22)
  
  Tissue1<-temp22$class[2,]
  Tissue2<-unique(temp22$class[2,])
  
  TF<-unique(temp22$class[3,])
  TF.n<-length(TF)
  
  temp3=file.path(output.file.dir,"Venn")
  
  if(!dir.exists(temp3))
  {
    dir.create(temp3)
  }
  
  for(i in 1:length(Tissue2)){
    po<-which(Tissue1 %in% Tissue2[i])
    
    print(po)
    
    if(length(po)==2)
    {
      png(file.path(temp3,paste0(paste0(colnames(temp22$class)[po],collapse = "-vs-"),".png")))
      dba.plotVenn(temp22,mask=po,main=Tissue2[i])
      dev.off()
      
    }else if(length(po)==4){
      po1<-po[c(1,3)]
      po2<-po[c(2,4)]
      
      png(file.path(temp3,paste0(paste0(colnames(temp22$class)[po1],collapse = "-vs-"),".png")))
      dba.plotVenn(temp22,mask=po1,main=Tissue2[i])
      dev.off()
      
      png(file.path(temp3,paste0(paste0(colnames(temp22$class)[po2],collapse = "-vs-"),".png")))
      dba.plotVenn(temp22,mask=po2,main=Tissue2[i])
      dev.off()
      
    }else
    {
      cat(paste0("For ",Tissue2[i],": Only one peak profile for ",TF.n," TFs\n"))
    }
  }

  p.common<-lapply(1:length(Tissue2),function(u,Tissue1,Tissue2,temp22){
    
    po<-which(Tissue1 %in% Tissue2[u])
    
    if(length(po)==2)
    {
      common.peaks<-dba.overlap(temp22,mask=po)
      y<-common.peaks$inAll
    }else if(length(po)==4){
      po1<-po[c(1,3)]
      po2<-po[c(2,4)]
      
      common.peaks.1<-dba.overlap(temp22,mask=po1)
      y1<-common.peaks.1$inAll
      
      common.peaks.2<-dba.overlap(temp22,mask=po2)
      y2<-common.peaks.2$inAll
      
      y<-list(y1=y1,y2=y2)
      
      names(y)[1]<-paste0(colnames(temp22$class)[po1],collapse = "-vs-")
      names(y)[2]<-paste0(colnames(temp22$class)[po2],collapse = "-vs-")
      
    }else{
      y<-NULL}
    y
  },Tissue1,Tissue2,temp22)
  
  names(p.common)<-Tissue2
  
  p.common<-unlist(p.common,recursive = F)
  
  p.common<-p.common[lapply(p.common,length) > 0]

  if(!dir.exists(file.path(output.file.dir,"common_peaks_bed")))
  {
    dir.create(file.path(output.file.dir,"common_peaks_bed"))
    dir.create(file.path(output.file.dir,"common_peaks_bed","ucsc"))
    dir.create(file.path(output.file.dir,"common_peaks_bed","igv"))
  }
  
  #output common peaks to bed files
  
  lapply(1:length(p.common),function(u,p.common,output.file.dir){
    x=p.common[[u]]

    x_name=names(p.common)[u]
    
    df <- data.frame(seqnames=seqnames(x),
                     #starts=start(x)-1,
                     starts=start(x),
                     ends=end(x),
                     names=c(rep(".", length(x))),
                     scores=elementMetadata(x)[,1],
                     strands=strand(x))

    #assign strand
    df.str <- data.frame(seqnames=seqnames(x),
                     #starts=start(x)-1,
                     starts=start(x),
                     ends=end(x),
                     names=c(rep(".", length(x))),
                     scores=elementMetadata(x)[,1],
                     strands=c(rep(".", length(x))))

    df.str.1<-df.str[-grep("random",df.str$seqnames),]
    
    df.str.2<-df.str.1

    df.str.3<-df.str.2[-grep("chrUn",df.str.2$seqnames),]

    write.table(df,file=file.path(output.file.dir,"common_peaks_bed",paste0(x_name,"_cp_with_header.bed")),
                col.names=TRUE,row.names = FALSE,quote=FALSE,sep="\t")

    write.table(df,file=file.path(output.file.dir,"common_peaks_bed","ucsc",paste0(x_name,"_4_ucsc.bed")),
                col.names=FALSE,row.names = FALSE,quote=FALSE,sep="\t")

    write.table(df,file=file.path(output.file.dir,"common_peaks_bed",paste0(x_name,"_common_peaks.bed")),
                col.names=FALSE,row.names = FALSE,quote=FALSE,sep="\t")

    write.table(df,file=file.path(output.file.dir,"common_peaks_bed","igv",paste0(x_name,"_4_igv.bed")),
                col.names=FALSE,row.names = FALSE,quote=FALSE,sep="\t")

    },p.common,output.file.dir)

  AnntationUsingChipSeeker(file.path(output.file.dir,"common_peaks_bed","igv"),"bed",file.path(output.file.dir,"common_peaks_bed")
                           ,txdb="hg19",DD=5000)
  
  return(p.common)
}