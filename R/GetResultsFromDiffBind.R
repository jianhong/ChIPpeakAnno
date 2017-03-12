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
    
    #pos=regexpr("\\.",x)
    
    #y=substr(x,1,pos-1)
    y=x
    y
  },sampID.v))
  
  #sampID.v.2<-substr(sampID.v,27,nchar(sampID.v)-3)
 
  colnames(temp$class)<-sampID.v.2
  #colnames(temp$class)<-sampID.v.2
  temp$class[1,]<-sampID.v.2
  
  #Merge the replicates of each set
  if(Mergereplicates=="yes"){
  temp2<-dba.peakset(temp,consensus = -DBA_REPLICATE)
  #temp22<-dba(temp2,mask =c(9,16,17:23))
  
  temp22<-dba(temp2,mask =c(17:24))
  }else{
    temp22<-temp
  }
  
  Tissue1<-temp22$class[2,]
  Tissue2<-unique(temp22$class[2,])
  
  TF<-unique(temp22$class[3,])
  TF.n<-length(TF)
  
  #Draw 
  #png(paste0(output.file.dir,"Venn.png"))
  
  #par(mfrow=c(2,3))
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
      
      #png(paste0(paste0(output.file.dir,"Venn/"),paste0(colnames(temp22$class)[po1],collapse = "-vs-"),".png"))
      
      png(file.path(temp3,paste0(paste0(colnames(temp22$class)[po1],collapse = "-vs-"),".png")))
      dba.plotVenn(temp22,mask=po1,main=Tissue2[i])
      dev.off()
      
      #png(paste0(paste0(output.file.dir,"Venn/"),paste0(colnames(temp22$class)[po2],collapse = "-vs-"),".png"))
      
      png(file.path(temp3,paste0(paste0(colnames(temp22$class)[po2],collapse = "-vs-"),".png")))
      dba.plotVenn(temp22,mask=po2,main=Tissue2[i])
      dev.off()
      
    }else
    {
      cat(paste0("For ",Tissue2[i],": Only one peak profile for ",TF.n," TFs\n"))
    }
  }
  #dev.off()
  #Get common peaks for TFs
  
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
  
  
  #library(EnsDb.Hsapiens.v75)
  #dd.hs<-toGRanges(EnsDb.Hsapiens.v75)


  #Tx <- transcripts(EnsDb.Hsapiens.v75)
  #dddd.hs<-Tx    # can not match

  # annotate each peak to genes

  # n.peaks<-length(temp$peaks)
  # lapply(1:n.peaks,function(u,temp,dd.hs,output.file.dir){
  #   x_name=colnames(temp$class)[u]
  #   x=toGRanges(temp$peaks[[u]])
  # 
  #   seqlevels(dd.hs,force=TRUE)<-seqinfo(x)@seqnames
  #   seqlevelsStyle(x) <- seqlevelsStyle(dd.hs)
  #   re.out.trimmed<-trim(x, use.names=TRUE)
  #   overlaps.anno<-annoPeaks(re.out.trimmed,dd.hs)
  # 
  #   write.table(overlaps.anno,file=paste0(output.file.dir,x_name,"_annotation_3_tr.txt"),row.names = FALSE,quote=FALSE,sep="\t")
  # },temp,dd.hs,output.file.dir)

  # if(!dir.exists(file.path(output.file.dir,"ucsc")))
  # {
  # dir.create(file.path(output.file.dir,"ucsc"))
  # }
  # 
  # if(!dir.exists(file.path(output.file.dir,"igv")))
  # {
  # dir.create(file.path(output.file.dir,"igv/"))
  # }
  
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
    
    #seqlevels(dd.hs,force=TRUE)<-seqinfo(x)@seqnames
    #seqlevelsStyle(x) <- seqlevelsStyle(dd.hs)
    #re.out.trimmed<-trim(x, use.names=TRUE)
    #overlaps.anno<-annoPeaks(re.out.trimmed,dd.hs)

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

    #df.str.2<-df.str.1[-which(as.character(df.str.1$starts)=="-1"),]

    df.str.2<-df.str.1

    df.str.3<-df.str.2[-grep("chrUn",df.str.2$seqnames),]

    write.table(df,file=file.path(output.file.dir,"common_peaks_bed",paste0(x_name,"_cp_with_header.bed")),col.names=TRUE,row.names = FALSE,quote=FALSE,sep="\t")

    write.table(df,file=file.path(output.file.dir,"common_peaks_bed","ucsc",paste0(x_name,"_4_ucsc.bed")),col.names=FALSE,row.names = FALSE,quote=FALSE,sep="\t")

    write.table(df,file=file.path(output.file.dir,"common_peaks_bed",paste0(x_name,"_common_peaks.bed")),col.names=FALSE,row.names = FALSE,quote=FALSE,sep="\t")

    write.table(df,file=file.path(output.file.dir,"common_peaks_bed","igv",paste0(x_name,"_4_igv.bed")),col.names=FALSE,row.names = FALSE,quote=FALSE,sep="\t")

    #write.table(x,file=paste0(output.file.dir,x_name,"_common_peaks.txt"),row.names = FALSE,quote=FALSE,sep="\t")

    #write.table(overlaps.anno,file=paste0(output.file.dir,x_name,"_annotation_2_tr.txt"),row.names = FALSE,quote=FALSE,sep="\t")
  },p.common,output.file.dir)


  AnntationUsingChipSeeker(file.path(output.file.dir,"common_peaks_bed","igv"),"bed",file.path(output.file.dir,"common_peaks_bed")
                           ,txdb="hg19",DD=5000)
  
  
  # overlap1.2<-dba.overlap(temp22,mask=c(1,2))
  # 
  # 
  # temp22$class[1,1]
  # temp22$class[2,2]
  # 
  # count.A=length(overlap1.2$onlyA@ranges@start)
  # count.B=length(overlap1.2$onlyB@ranges@start)
  # count.All=length(overlap1.2$inAll@ranges@start)
  # 
  # cat(count.A," ",count.B," ",count.All,"\n")
  # 
  # .Weight<-c(0,count.A,count.B,count.All)
  # A=c(0,1,0,1)
  # B=c(0,0,1,1)
  # tt<-as.matrix(cbind(A,B,.Weight))
  # 
  # row.names(tt)=c("00","10","01","11")
  # colnames(tt)[1:2]=c(temp22$class[1,1],temp22$class[1,2])
  # 
  # Vstem.test <- Venn()
  # Vstem.test@IndicatorWeight<-tt
  # Cstem3.test <- compute.Venn(Vstem.test,doWeights=TRUE)
  # c.3.set=sum(Cstem3.test@IndicatorWeight[,3])
  # Cstem3.test@IndicatorWeight[which(row.names(Cstem3.test@IndicatorWeight)=="00"),3]=6000-c.3.set
  # 
  # #SetLabels <- VennGetSetLabels(Cstem3.test)
  # #FaceLabels <- VennGetFaceLabels(Cstem3.test)
  # 
  # UseFisher <- function(temp.ct,index.A,index.B,totalN) {
  #   total.peaks=totalN
  #   A=sum(temp.ct[which(temp.ct[,index.A]==1&temp.ct[,index.B]==1),3])
  #   B=sum(temp.ct[which(temp.ct[,index.A]==1&temp.ct[,index.B]==0),3])
  #   C=sum(temp.ct[which(temp.ct[,index.A]==0&temp.ct[,index.B]==1),3])
  #   D=total.peaks-(A+B+C)
  #   ctb<-matrix(c(A,B,C,D),nrow = 2,dimnames =list(c("In", "Not"),c("In", "Not")))
  # 
  #   #re<-fisher.test(ctb)
  #   print(ctb)
  #   re.fisher<-fisher.test(ctb, alternative='greater')[c("p.value","estimate")]
  #   re.fisher
  # }
  # 
  # index.A<-grep(temp22$class[1,1],colnames(tt))
  # index.B<-grep(temp22$class[1,2],colnames(tt))
  # tempRe<-UseFisher(tt,index.A,index.B,6000)
  # 
  # cat(tempRe$p.value," ",tempRe$estimate)
  # 
  # #Visualization of peak sets overlapping using Area-Proportional Venn diagrams
  # 
  # png(paste0(output.file.dir,"APVenn.png"))
  # 
  # plot(Cstem3.test)
  # x <- 0.35
  # y <- 0.15
  # rot <- 0
  # grid.text(paste0("fisher.test p value: ",tempRe$p.value), x=x, y=y, rot=rot,gp=gpar(fontsize=15, col="blue"))
  # 
  # dev.off()
  
  #Vstem.test@IndicatorWeight<-Vstem@IndicatorWeight
  #Vstem.test@IntersectionSets<-Vstem.test@IntersectionSets
  # temp3<-dba.peakset(temp2,bRetrieve = TRUE)
  # temp.no.replicate<-dba(dba.peakset(temp2,c(9,16),bRetrieve=TRUE,DataTyp=DBA_DATA_FRAME))
  # temp.replicate<-dba.peakset(temp2,17:23,bRetrieve=TRUE,DataTyp=DBA_DATA_FRAME)
  # temp.c<-dba.peakset(peaks = temp.replicate)
  # temp3<-dba.peakset(temp2,temp2$masks$Consensus,bRetrieve = TRUE)
  # temp3<-dba.peakset(NULL,peaks=temp2$masks$Consensus)
  # 
  
  return(p.common)
}