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


  test <- dba.peakset(NULL,peaks="/media/H_driver/2016/Yang/Test/A.bed",
                      peak.caller="bed", sampID="A",tissue="brain",
                      factor="ER",condition="A",replicate=1)

  test <- dba.peakset(test,peaks="/media/H_driver/2016/Yang/Test/B.bed",
                      peak.caller="bed", sampID="B",tissue="brain",
                      factor="ER",condition="B",replicate=1)

  test<-dba.peakset(test,consensus=DBA_FACTOR, minOverlap=2)

  test.OL <- dba.overlap(test,test$masks$Consensus)

  dba.plotVenn(test,c(1,2))


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

  ctb<-matrix(c(13289,1244, 7899,0),nrow = 2,dimnames =list(c("InKO_BM_SMC1", "NotKO_BM_SMC1"),c("In_KO_BM_Rad21", "Not_KO_BM_Rad21")))

  #ctb.2<-matrix(c(13289,7899,13289,1244),nrow = 2,dimnames =list(c("Overlaped", "Not overlaped"),c("KO_BM_SMC1", "KO_BM_Rad21")))

  re<-fisher.test(ctb)

  re.2<-fisher.test(ctb.2)


  WT.triple <- dba.peakset(NULL,peaks=file.name.3[[2]],
                      peak.caller="bed", sampID=names(file.name.3)[2],tissue="Yang",
                      factor="ER",condition="WT",replicate=1)

  WT.triple <- dba.peakset(WT.triple,peaks=file.name.3[[3]],
                      peak.caller="bed", sampID=names(file.name.3)[3],tissue="Yang",
                      factor="ER",condition="WT",replicate=1)

  WT.triple <- dba.peakset(WT.triple,peaks=file.name.3[[5]],
                      peak.caller="bed", sampID=names(file.name.3)[5],tissue="Yang",
                      factor="ER",condition="WT",replicate=1)

  Binding.ASXL1<-WT.triple$binding[which(WT.triple$binding[,4]==1),]

  WT.triple.binding<-WT.triple$binding[which(WT.triple$binding[,4]==1&WT.triple$binding[,5]==1&WT.triple$binding[,6]==1),]

  dba.plotVenn(WT.triple,c(2,3,1))

  WT.triple.binding[,1]<-paste0("chr",WT.triple.binding[,1])

  write.table(WT.triple.binding[,1:3],file=paste0(out.dir.name,"WT_triple_overlap.bed"),sep="\t",
              quote = FALSE,row.names = FALSE,col.names = FALSE)

  # sample.name<-sapply(strsplit(names(re.out),split="_peaks_"),"[[",1)
  ##
  # names(re.out)=sample.name
  #
  # A=re.out[[1]] #10_2470IUPUI_WT_BM_ASXL1
  # B=re.out[[2]] #10_WT_BM_ASXL1
  # C=re.out[[3]] #11_2470IUPUI_WT_BM_SMC1
  # D=re.out[[5]] #13_2470IUPUI_WT_BM_Rad21
  #
  # A.name=sample.name[1] #10_2470IUPUI_WT_BM_ASXL1
  # B.name=sample.name[2] #10_WT_BM_ASXL1
  # C.name=sample.name[3] #11_2470IUPUI_WT_BM_SMC1
  # D.name=sample.name[5] #13_2470IUPUI_WT_BM_Rad21

  #A <- createRandomRegions(nregions=50, length.mean=5000000, length.sd=3000000)
  #B <- c(A[1:25], createRandomRegions(nregions=25, length.mean=500000, length.sd=30000))
  #class(A)
  #class(B)


  # TestOverlap <- function(A, B, A.name,B.name,out.dir.name,n_permutations) {
  #   NumOverLap<-numOverlaps(A, B, count.once=TRUE)
  #   pt <- overlapPermTest(A=A, B=B, ntimes=n_permutations)
  #   png(filename=paste0(out.dir.name,A.name,B.name,n_permutations,"_overlap_test.png"))
  #   plot(pt)
  #   dev.off()
  #
  #   lz <- localZScore(pt=pt, A=A, B=B)
  #   png(filename=paste0(out.dir.name,A.name,B.name,n_permutations,"_localZScore.png"))
  #   plot(lz)
  #   dev.off()
  # }
  #
  # TestOverlap(A,B,A.name,B.name,out.dir.name,100)
  # TestOverlap(A,C,A.name,C.name,out.dir.name,100)
  # TestOverlap(A,D,A.name,D.name,out.dir.name,100)
  #
  # TestOverlap(B,C,B.name,C.name,out.dir.name,100)
  # TestOverlap(B,D,B.name,D.name,out.dir.name,100)
  # TestOverlap(C,D,C.name,D.name,out.dir.name,100)
  #

 # temp.name<-strsplit(names(file.name.2),split="\\.")

 #a = readBed(system.file("extdata", "examples/combined_regions.bed",
 #                        package="LOLA"))


  # venn.plot <- venn.diagram(
  #   x = re.out[c(1,2)],
  #   filename = paste0(out.dir.name,names(re.out)[1],names(re.out)[2],"_overlap_venn.tiff"),
  #   height = 3000,
  #   width = 3500,
  #   resolution = 1000,
  #   col = "black",
  #   lty = "dotted",
  #   lwd = 1,
  #   fill = c("red","blue"),
  #   alpha = 0.50,
  #   label.col = c(rep("white",3)),
  #   cex = 0.5,
  #   fontfamily = "serif",
  #   fontface = "bold",
  #   cat.col = c("red","blue"),
  #   cat.cex = 0.5,
  #   cat.pos = 0.5,
  #   cat.dist = 0.05,
  #   cat.fontfamily = "serif"
  # )
  #
  # venn.plot <- venn.diagram(
  #   x = re.out[c(1,3,5)],
  #   filename =paste0(out.dir.name,names(re.out)[1],names(re.out)[3],names(re.out)[5],"_overlap_venn.tiff"),
  #   col = "black",
  #   lty = "dotted",
  #   lwd = 2,
  #   fill = c("red", "orange", "blue"),
  #   alpha = 0.50,
  #   label.col = c(rep("white",7)),
  #   cex = 1,
  #   fontfamily = "serif",
  #   fontface = "bold",
  #   cat.col = c("red", "orange", "blue"),
  #   cat.cex = 0.8,
  #   cat.fontfamily = "serif"
  # )


  return(re.out)
}
