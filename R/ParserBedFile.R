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
#' out.dir.name="/media/H_driver/2016/Yang/Results/"
#' output.test.file= "SMC1A_RAD21_overlap_test"
#'
#' ParserBedFile(dir.name,input.file.pattern,out.dir.name,output.test.file)
#'
ParserBedFile<-function(dir.name,input.file.pattern,out.dir.name,output.test.file){

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

  A=re.out[[2]]
  B=re.out[[4]]

  #A <- createRandomRegions(nregions=50, length.mean=5000000, length.sd=3000000)
  #B <- c(A[1:25], createRandomRegions(nregions=25, length.mean=500000, length.sd=30000))
  #class(A)
  #class(B)

  NumOverLap<-numOverlaps(A, B, count.once=TRUE)
  pt <- overlapPermTest(A=A, B=B, ntimes=50)

  png(filename=paste0(out.dir.name,output.test.file,".png"))
  plot(pt)
  dev.off()

 # temp.name<-strsplit(names(file.name.2),split="\\.")

 #a = readBed(system.file("extdata", "examples/combined_regions.bed",
 #                        package="LOLA"))
}
