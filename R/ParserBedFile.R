#' Title
#'
#' @param dir.name
#' @param input.file.pattern
#'
#' @return
#' @export
#'
#' @examples
#'
#' dir.name="/media/H_driver/2016/Yang/BedFromPeakCall/"
#' input.file.pattern="*bam_mm_summits_from_PeakCall4Yang.bed"
#'
#' out.dir.name="/media/H_driver/2016/Yang/Results/"
#' output.test.file= "SMC1A_RAD21_overlap_test"
#'
#'
#' dir.name.2="/media/H_driver/2016/Yang/BedFromPeakCall_mm_100e_07/"
#' input.file.pattern.2="*_bam_mm_1_from_PeakCall4Yang.bed"
#'
#' out.dir.name.2="/media/H_driver/2016/Yang/Results/"
#' output.test.file.2= "SMC1A_RAD21_mm_100e-07_peaks_overlap_test"
#'
#' ParserBedFile(dir.name,input.file.pattern,out.dir.name)
#'
#' re.mm.100e.07<-ParserBedFile(dir.name.2,input.file.pattern.2,out.dir.name)
#'
#'
#' save.image(file=paste0(out.dir.name,"re_save_2.RData"))
#' savehistory(file=paste0(out.dir.name,"re_save_2.Rhistory"))
#'
ParserBedFile<-function(dir.name,input.file.pattern,out.dir.name){

  file.name=paste0(dir.name,dir(dir.name,recursive = TRUE,pattern=input.file.pattern))
  file.name.2<-as.list(file.name)

  names(file.name.2)=sapply(strsplit(file.name,split="\\/"),"[[",7)

  print(file.name.2)

  file.name.3<-file.name.2[-6]

  re.out<-lapply(file.name.3,function(u){
    re=readBed(u)
    #colnames(re)=c("Count","GeneName")
    re
  })

  sample.name<-sapply(strsplit(sapply(strsplit(names(re.out),split="_i"),"[[",1),split="2470IUPUI_"),"[[",2)

  A=re.out[[2]] #SMC1
  B=re.out[[4]] #Rad21
  C=re.out[[1]] #ASXL1

  A.name=sample.name[2] #SMC1 name
  B.name=sample.name[4] #Rad21 name
  C.name=sample.name[1] #ASXL1 name

  #A <- createRandomRegions(nregions=50, length.mean=5000000, length.sd=3000000)
  #B <- c(A[1:25], createRandomRegions(nregions=25, length.mean=500000, length.sd=30000))
  #class(A)
  #class(B)


  TestOverlap <- function(A, B, A.name,B.name,out.dir.name) {
    NumOverLap<-numOverlaps(A, B, count.once=TRUE)
    pt <- overlapPermTest(A=A, B=B, ntimes=50)
    png(filename=paste0(out.dir.name,A.name,B.name,"_overlap_test.png"))
    plot(pt)
    dev.off()

    lz <- localZScore(pt=pt, A=A, B=B)
    png(filename=paste0(out.dir.name,A.name,B.name,"_localZScore.png"))
    plot(lz)
    dev.off()
  }

  TestOverlap(A,B,A.name,B.name,out.dir.name)
  TestOverlap(A,C,A.name,C.name,out.dir.name)
  TestOverlap(B,C,B.name,C.name,out.dir.name)

 # temp.name<-strsplit(names(file.name.2),split="\\.")

 #a = readBed(system.file("extdata", "examples/combined_regions.bed",
 #                        package="LOLA"))

  return(re.out)
}
