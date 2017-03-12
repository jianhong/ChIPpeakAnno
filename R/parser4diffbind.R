#' parser4diffbind
#' 
#' parser4diffbind
#'
#' @param input.sample.file
#' @param input.bed.dir
#' @param input.file.pattern
#'
#' @return
#' @export
#'
#' @examples
#'  
#' input.sample.file = "/Volumes/Bioinformatics$/2017/DannyNewData/SampleID_INFO_ChIP_new_Danny.csv"
#' input.bed.dir = "/Volumes/Bioinformatics$/2017/DannyNewData/PeakCall/"
#' input.file.pattern = "bed"
#' resmcf<-parser4diffbind(input.sample.file,input.bed.dir,input.file.pattern)
#'
parser4diffbind<-function(input.sample.file,input.bed.dir,input.file.pattern){

  
  sample.info<-read.csv(input.sample.file,header = TRUE)
  sample.info.2 <- sample.info[,-dim(sample.info)[2]]
  
  out.dir.name = dirname(input.sample.file)
  
  re<-ParserReadFiles(input.bed.dir,input.file.pattern,out.dir.name)
  
  re.bed<-re$input
  
  re.peaks.only.bed.2<-FL(re.bed,'peaks')
  re.summits.only.bed<-FL(re.bed,'summits')
  
  
  #put peak files into database
  re.peaks.only.bed.3<-list_to_df(re.peaks.only.bed.2)
  
  pos <- regexpr('R1', re.peaks.only.bed.3$ID)
  ID2 <- substr(re.peaks.only.bed.3$ID,1,pos+1)
  
  re.peaks.only.bed.3 <- cbind(ID2,re.peaks.only.bed.3)
  
  colnames(re.peaks.only.bed.3)=c("ID","ID2") 
  
  re.peaks.only.bed.4<-merge(sample.info.2,re.peaks.only.bed.3,by="ID")

  replicate<-apply(re.peaks.only.bed.4,1,function(x){
     ID=as.character(x[1])
     cell=as.character(x[2])
     TF=as.character(x[3])
     ID2=as.character(x[4])
     p1=nchar(ID)+2
     p2=regexpr(TF,ID2)
     p2=p2-2
     xx=substr(ID2,p1,p2)
     cell=substr(xx,1,nchar(xx)-2)
     rep=substr(xx,nchar(xx),nchar(xx))
     file=as.character(x[6])
     condition=paste(cell,TF,sep=":")
     z <- data.frame(t(c(ID2,cell,TF,condition,rep,file)))
     z
   })
   
   res <- do.call(rbind.data.frame,replicate)
   colnames(res) <- c("ID","Cell","Ab","Cell_Ab","Replicate","File")
   

  re.peaks.only.bed.5 <- res
  
  mcf<-NULL

  for(i in 1:dim(re.peaks.only.bed.5)[1]){

    peaks.v=readBed(as.character(re.peaks.only.bed.5[i,6]))

    peak.caller.v="macs142"
    sampID.v=as.character(re.peaks.only.bed.5[i,1])

    sampID.v=substr(sampID.v,1,nchar(sampID.v)-5)

    tissue.v=as.character(re.peaks.only.bed.5[i,2])
    factor.v=as.character(re.peaks.only.bed.5[i,3])
    condition.v=as.character(re.peaks.only.bed.5[i,4])
    
    replicate.v=as.numeric(as.character(re.peaks.only.bed.5[i,5]))

    if(is.null(mcf)){
      mcf <- dba.peakset(NULL,peaks=peaks.v,peak.caller=peak.caller.v, sampID=sampID.v,tissue=tissue.v,
                          factor=factor.v,condition=condition.v,replicate=replicate.v)
    }else{
      mcf <- dba.peakset(mcf,peaks=peaks.v,peak.caller=peak.caller.v, sampID=sampID.v,tissue=tissue.v,
                          factor=factor.v,condition=condition.v,replicate=replicate.v)
    }
  }

  print(mcf)
  
  temp3=file.path(re$output,"BindDiff")
  if(!dir.exists(temp3)){dir.create(temp3)}
   
  save(mcf,file=file.path(temp3,"mcf.RData"))
  
   pdf(file.path(temp3,"corrmap.pdf"))
   dba.plotHeatmap(mcf,margin=23)
   dev.off()
 
   GetResultsFromDiffBind(mcf,"yes",temp3)
  
   return(mcf)
}
