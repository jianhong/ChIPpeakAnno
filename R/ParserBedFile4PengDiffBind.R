#' ParserBedFile4PengDiffBind
#'
#' @param dir.name
#' @param input.file.pattern
#'
#' @return
#' @export
#'
#' @examples
#'
#' dir.name="/media/H_driver/2016/Yang/BedFromPeng/"
#' input.file.pattern="*peaks_from_PeakCall4Yang.bed"
#' out.dir.name="/media/H_driver/2016/Yang/Results/"

#'
#' re.from.bed.peng.4.venn<-ParserBedFile4PengVenn(dir.name,input.file.pattern,out.dir.name)
#'
#'
#' save.image(file=paste0(out.dir.name,"re_save_peng.RData"))
#' savehistory(file=paste0(out.dir.name,"re_save_peng.Rhistory"))
#'

ParserBedFile4PengDiffBind<-function(dir.name,input.file.pattern,out.dir.name){

  file.name=paste0(dir.name,dir(dir.name,recursive = TRUE,pattern=input.file.pattern))
  file.name.2<-as.list(file.name)

  names(file.name.2)=sapply(strsplit(file.name,split="\\/"),"[[",7)

  print(file.name.2)

  file.name.3<-file.name.2

  sample.name<-sapply(strsplit(names(file.name.3),split="_peaks_"),"[[",1)

  names(file.name.3)=sample.name


  # test <- dba.peakset(NULL,peaks="/media/H_driver/2016/Yang/Test/A.bed",
  #                     peak.caller="bed", sampID="A",tissue="brain",
  #                     factor="ER",condition="A",replicate=1)
  # 
  # test <- dba.peakset(test,peaks="/media/H_driver/2016/Yang/Test/B.bed",
  #                     peak.caller="bed", sampID="B",tissue="brain",
  #                     factor="ER",condition="B",replicate=1)
  # 
  # test<-dba.peakset(test,consensus=DBA_FACTOR, minOverlap=2)
  # 
  # test.OL <- dba.overlap(test,test$masks$Consensus)
  # 
  # dba.plotVenn(test,c(1,2))


  mcf7 <- dba.peakset(NULL,peaks=file.name.3[[1]],
                    peak.caller="bed", sampID=names(file.name.3)[1],tissue="MCF7",
                    factor="ER",condition="WT",replicate=1)

  mcf7 <- dba.peakset(mcf7,peaks=file.name.3[[2]],
                      peak.caller="bed", sampID=names(file.name.3)[2],tissue="MCF7",
                      factor="ER",condition="WT",replicate=1)


  mcf7 <- dba.peakset(mcf7,peaks=file.name.3[[3]],
                      peak.caller="bed", sampID=names(file.name.3)[3],tissue="MCF7",
                      factor="ER",condition="WT",replicate=1)

  mcf7 <- dba.peakset(mcf7,peaks=file.name.3[[4]],
                      peak.caller="bed", sampID=names(file.name.3)[4],tissue="MCF7",
                      factor="ER",condition="KO",replicate=1)

  mcf7 <- dba.peakset(mcf7,peaks=file.name.3[[5]],
                      peak.caller="bed", sampID=names(file.name.3)[5],tissue="MCF7",
                      factor="ER",condition="WT",replicate=1)

  mcf7 <- dba.peakset(mcf7,peaks=file.name.3[[6]],
                      peak.caller="bed", sampID=names(file.name.3)[6],tissue="MCF7",
                     factor="ER",condition="KO",replicate=1)

  #Remake peak dba
  Peng <- dba.peakset(NULL,peaks=file.name.3[[2]],
                      peak.caller="bed", sampID=names(file.name.3)[2],tissue="Unkown",
                      factor="ASXL1",condition="WT",replicate=1)

  Peng <- dba.peakset(Peng,peaks=file.name.3[[3]],
                      peak.caller="bed", sampID=names(file.name.3)[3],tissue="Unkown",
                      factor="SMC1",condition="WT",replicate=1)

  Peng <- dba.peakset(Peng,peaks=file.name.3[[4]],
                      peak.caller="bed", sampID=names(file.name.3)[4],tissue="Unkown",
                      factor="SMC1",condition="KO",replicate=1)

  Peng <- dba.peakset(Peng,peaks=file.name.3[[5]],
                      peak.caller="bed", sampID=names(file.name.3)[5],tissue="Unkown",
                      factor="Rad21",condition="WT",replicate=1)

  Peng <- dba.peakset(Peng,peaks=file.name.3[[6]],
                      peak.caller="bed", sampID=names(file.name.3)[6],tissue="Unkown",
                      factor="Rad21",condition="KO",replicate=1)


  Peng

  olap.rate <- dba.overlap(mcf7,mode=DBA_OLAP_RATE)

  olap.rate

  dba.overlap(mcf7,mode=DBA_OLAP_RATE)
  #dba.plotVenn(mcf7,c(1,3,5))

  #dba.plotVenn(mcf7,c(2,3,5))

  mcf7.consensus.3.5 <- dba.peakset(mcf7,c(3,5), minOverlap=2, bRetrieve=TRUE)

  mcf7<- dba.peakset(mcf7,c(3,5), minOverlap=2, sampID = "WT_SMC1A_RAD21")

  mcf7<- dba.peakset(mcf7,c(4,6), minOverlap=2, sampID = "KO_SMC1A_RAD21")

  dba.plotVenn(mcf7,c(9,10))


  dba.plotVenn(mcf7,c(3,5,2))

  dba.plotVenn(mcf7,c(3,5))

  dba.plotVenn(mcf7,c(4,6))

  #ctb<-matrix(c(13289,1244, 7899,0),nrow = 2,dimnames =list(c("InKO_BM_SMC1", "NotKO_BM_SMC1"),c("In_KO_BM_Rad21", "Not_KO_BM_Rad21")))

  universe.count=
  d=
  ctb2<-matrix(c(13289,1244,7899,d),nrow = 2,dimnames =list(c("InKO_BM_SMC1", "NotKO_BM_SMC1"),c("In_KO_BM_Rad21", "Not_KO_BM_Rad21")))

  #re<-fisher.test(ctb)
  ctb2
  fisher.test(ctb2, alternative='greater')[c("p.value","estimate")]


  a=2309
  b=21991-2309
  c=2762-2309
  d=32381-(a+b+c)

  ctb3<-matrix(c(a,b,c,d),nrow = 2,dimnames =list(c("InKO_BM_SMC1", "NotKO_BM_SMC1"),c("In_KO_BM_Rad21", "Not_KO_BM_Rad21")))

  #re<-fisher.test(ctb)
  ctb3
  fisher.test(ctb2, alternative='greater')[c("p.value","estimate")]


  #Make peak dab for WT only
  WT.triple <- dba.peakset(NULL,peaks=file.name.3[[2]],
                      peak.caller="bed", sampID=names(file.name.3)[2],tissue="Unknown",
                      factor="ASXL1",condition="WT",replicate=1)

  WT.triple <- dba.peakset(WT.triple,peaks=file.name.3[[3]],
                      peak.caller="bed", sampID=names(file.name.3)[3],tissue="Unknown",
                      factor="SMC1",condition="WT",replicate=1)

  WT.triple <- dba.peakset(WT.triple,peaks=file.name.3[[5]],
                      peak.caller="bed", sampID=names(file.name.3)[5],tissue="Unknown",
                      factor="Rad21",condition="WT",replicate=1)

  dba.plotVenn(WT.triple,c(2,3,1))


  WT.SMC1.RAD21 <- dba.peakset(NULL,peaks=file.name.3[[3]],
                           peak.caller="bed", sampID=names(file.name.3)[3],tissue="Unknown",
                           factor="SMC1",condition="WT",replicate=1)

  WT.SMC1.RAD21 <- dba.peakset(WT.SMC1.RAD21,peaks=file.name.3[[5]],
                           peak.caller="bed", sampID=names(file.name.3)[5],tissue="Unknown",
                           factor="Rad21",condition="WT",replicate=1)

  WT.SMC1.RAD21.binding<-WT.SMC1.RAD21$binding[,1:3]

  WT.SMC1.RAD21.binding[,1]<-paste0("chr",WT.SMC1.RAD21.binding[,1])

  WT.SMC1.RAD21.binding.dba<-dba.peakset(NULL,peaks= paste0(out.dir.name,"SMC1_RAD21_overlap.bed"),
              peak.caller="bed", sampID="SMC1_RAD21",tissue="Unknown",
              factor="SMC1_RAD21",condition="WT",replicate=1)

  WT.SMC1.RAD21.binding[,1]<-paste0("chr",WT.SMC1.RAD21.binding[,1])

  write.table(WT.SMC1.RAD21.binding,file=paste0(out.dir.name,"SMC1_RAD21_overlap.bed"),sep="\t",
              quote = FALSE,row.names = FALSE,col.names = FALSE)

  WT.SMC1.RAD21.binding.dba<-dba.peakset(NULL,peaks= paste0(out.dir.name,"SMC1_RAD21_overlap.bed"),
                                         peak.caller="bed", sampID="SMC1_RAD21",tissue="Unknown",
                                         factor="SMC1_RAD21",condition="WT",replicate=1)

  WT.SMC1.RAD21.binding.dba<-dba.peakset(WT.SMC1.RAD21.binding.dba,peaks= file.name.3[[2]],
                                     peak.caller="bed", sampID=names(file.name.3)[2],tissue="Unknown",
                                     factor="ASXL1",condition="WT",replicate=1)


  WT.SMC1.RAD21.binding.dba

  dba.plotVenn(WT.SMC1.RAD21.binding.dba,c(1,2))

  WT.SMC1.Yes<-WT.triple$binding[which(WT.triple$binding[,4]==1),]
  WT.SMC1.No<-WT.triple$binding[which(WT.triple$binding[,4]==0),]

  dim(WT.SMC1.Yes)

  FreqFind <- function(WT.SMC1.Yes) {
    A=dim(WT.SMC1.Yes[which(WT.SMC1.Yes[,5]==1&WT.SMC1.Yes[,6]==1),])[1]
    B=dim(WT.SMC1.Yes[which(WT.SMC1.Yes[,5]==1&WT.SMC1.Yes[,6]==0),])[1]
    C=dim(WT.SMC1.Yes[which(WT.SMC1.Yes[,5]==0&WT.SMC1.Yes[,6]==1),])[1]
    D=22159-dim(WT.SMC1.Yes[which(WT.SMC1.Yes[,5]==0&WT.SMC1.Yes[,6]==0),])[1]
    re<-c(A,B,C,D)
    re
  }

  AA=FreqFind(WT.SMC1.Yes)
  BB=FreqFind(WT.SMC1.No)

  mytable <- array(c(AA,
                     BB),
                   dim = c(2,2,2),
                   dimnames = list(
                     Is_C = c('Yes','No'),
                     Is_B = c('Yes','No'),
                     Is_A = c('Yes','No')))

  mytable
  mantelhaen.test(mytable,exact = TRUE, alternative = "greater")


  dba.plotVenn(WT.triple,c(1,2,3))

  WT.triple.binding<-WT.triple$binding[which(WT.triple$binding[,4]==1&WT.triple$binding[,5]==1&WT.triple$binding[,6]==1),]


  Binding.ASXL1<-WT.triple$binding[which(WT.triple$binding[,4]==1),]

  WT.triple.binding<-WT.triple$binding[which(WT.triple$binding[,4]==1&WT.triple$binding[,5]==1&WT.triple$binding[,6]==1),]

  dba.plotVenn(WT.triple)

  WT.triple.binding[,1]<-paste0("chr",WT.triple.binding[,1])

  write.table(WT.triple.binding[,1:3],file=paste0(out.dir.name,"WT_triple_overlap.bed"),sep="\t",
              quote = FALSE,row.names = FALSE,col.names = FALSE)


  return(re.out)
}
