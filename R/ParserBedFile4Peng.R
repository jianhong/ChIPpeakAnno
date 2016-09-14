#' ParserBedFile4Peng
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
#' input.file.pattern="*.bed"
#'
#' out.dir.name="/media/H_driver/2016/Yang/Results/"

#'
#' re.from.bed.peng<-ParserBedFile4Peng(dir.name,input.file.pattern,out.dir.name)
#'
#'
#' save.image(file=paste0(out.dir.name,"re_save_peng.RData"))
#' savehistory(file=paste0(out.dir.name,"re_save_peng.Rhistory"))
#'
ParserBedFile4Peng<-function(dir.name,input.file.pattern,out.dir.name){

  file.name=paste0(dir.name,dir(dir.name,recursive = TRUE,pattern=input.file.pattern))
  file.name.2<-as.list(file.name)

  names(file.name.2)=sapply(strsplit(file.name,split="\\/"),"[[",7)

  print(file.name.2)

  #file.name.3<-file.name.2[-6]

  re.out<-lapply(file.name.2,function(u){
    re=readBed(u)
    #colnames(re)=c("Count","GeneName")
    re
  })

  sample.name<-sapply(strsplit(names(re.out),split="_peaks_"),"[[",1)

  A=re.out[[1]] #10_2470IUPUI_WT_BM_ASXL1
  B=re.out[[2]] #10_WT_BM_ASXL1
  C=re.out[[3]] #11_2470IUPUI_WT_BM_SMC1
  D=re.out[[5]] #13_2470IUPUI_WT_BM_Rad21

  A.name=sample.name[1] #10_2470IUPUI_WT_BM_ASXL1
  B.name=sample.name[2] #10_WT_BM_ASXL1
  C.name=sample.name[3] #11_2470IUPUI_WT_BM_SMC1
  D.name=sample.name[5] #13_2470IUPUI_WT_BM_Rad21

  #A <- createRandomRegions(nregions=50, length.mean=5000000, length.sd=3000000)
  #B <- c(A[1:25], createRandomRegions(nregions=25, length.mean=500000, length.sd=30000))
  #class(A)
  #class(B)


  TestOverlap <- function(A, B, A.name,B.name,out.dir.name,n_permutations) {
    NumOverLap<-numOverlaps(A, B, count.once=TRUE)
    pt <- overlapPermTest(A=A, B=B, ntimes=n_permutations)
    png(filename=paste0(out.dir.name,A.name,B.name,n_permutations,"_overlap_test.png"))
    plot(pt)
    dev.off()

    lz <- localZScore(pt=pt, A=A, B=B)
    png(filename=paste0(out.dir.name,A.name,B.name,n_permutations,"_localZScore.png"))
    plot(lz)
    dev.off()
  }

  TestOverlap(A,B,A.name,B.name,out.dir.name,100)
  TestOverlap(A,C,A.name,C.name,out.dir.name,100)
  TestOverlap(A,D,A.name,D.name,out.dir.name,100)

  TestOverlap(B,C,B.name,C.name,out.dir.name,100)
  TestOverlap(B,D,B.name,D.name,out.dir.name,100)
  TestOverlap(C,D,C.name,D.name,out.dir.name,100)


 # temp.name<-strsplit(names(file.name.2),split="\\.")

 #a = readBed(system.file("extdata", "examples/combined_regions.bed",
 #                        package="LOLA"))

  return(re.out)
}
