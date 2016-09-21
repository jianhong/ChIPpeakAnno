#' ParserBedFile4PengXls
#'
#' @param dir.name
#' @param input.file.pattern
#'
#' @return
#' @export
#'
#' @examples
#'
#' dir.name="/media/H_driver/2016/Yang/MACS/MACS/"
#' input.file.pattern="*.xls"
#'
#' out.dir.name="/media/H_driver/2016/Yang/Results/"

#'
#' re.from.bed.peng.Xls<-ParserBedFile4PengXls(dir.name,input.file.pattern,out.dir.name)
#'
#' re.from.bed.peng.Xls.2<-ParserBedFile4PengXls(dir.name,input.file.pattern,out.dir.name)
#'
#' re.from.bed.peng.Xls.3<-ParserBedFile4PengXls(dir.name,input.file.pattern,out.dir.name)
#'
#' save.image(file=paste0(out.dir.name,"re_save_peng.RData"))
#' savehistory(file=paste0(out.dir.name,"re_save_peng.Rhistory"))
#'
ParserBedFile4PengXls<-function(dir.name,input.file.pattern,out.dir.name){

  dir.name=reformatPath(dir.name)
  out.dir.name=reformatPath(out.dir.name)
  
  file.name=paste0(dir.name,dir(dir.name,recursive = TRUE,pattern=input.file.pattern))
  file.name.2<-as.list(file.name)

  names(file.name.2)=sapply(strsplit(file.name,split="\\/"),"[[",8)

  print(file.name.2)

  re.out<-lapply(names(file.name.2),function(u){


    #print(names(u))
    uu=file.name.2[[u]]

    re=read.table(uu,header=TRUE)
    colnames(re)[c(7,9)]=c("-10*LOG10(pvalue)","FDR(%)")
    temp<-re
    temp$start<-temp$start-1
    temp$end<-temp$end-1
    temp$summit<-temp$summit-1
    temp

    summitPeak<-temp$start+temp$summit

    temp2<-temp

    temp2$start<-summitPeak-49
    temp2$end<-summitPeak+50
    temp2$summit<-summitPeak
    temp2$length<-temp2$end-temp2$start+1
    temp2
    rownames(temp2)<-paste0("MACS_peak_",rownames(temp2))
    temp3<-toGRanges(temp2)
    names(temp3)<-rownames(temp2)


    dd.GRCm39.mm10<-toGRanges(EnsDb.Mmusculus.v75)
    genome(temp3)<-genome(dd.GRCm39.mm10)

    seqlevels(dd.GRCm39.mm10,force=TRUE) <- seqlevels(temp3)

    seqinfo(temp3)<-seqinfo(dd.GRCm39.mm10)

    temp3.trimmed<-trim(temp3, use.names=TRUE)

    genome.mm10<-getBSgenome("BSgenome.Mmusculus.UCSC.mm10")

    genome(temp3.trimmed)<-"mm10"

    seq.temp3<-getAllPeakSequence(temp3.trimmed,genome=genome.mm10)

    write2FASTA(seq.temp3, paste0(out.dir.name,u,".fa"))

    re2<-list(originalPeak=re,down1Peak=temp,aroundSummit100Peak=temp2,GR=temp3)
    re2
    })

  sample.name<-sapply(strsplit(names(file.name.2),split="_peaks_"),"[[",1)

  names(re.out)<-sample.name


  return(re.out)
}
