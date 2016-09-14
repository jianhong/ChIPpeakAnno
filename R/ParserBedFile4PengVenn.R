#' ParserBedFile4PengVenn
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
#' input.file.pattern="*combine_chr_st_end.bed"
#'
#' input.file.pattern="*combine_chr_st.bed"
#'
#'
#'
#' out.dir.name="/media/H_driver/2016/Yang/Results/"

#'
#' re.from.bed.peng.4.venn<-ParserBedFile4PengVenn(dir.name,input.file.pattern,out.dir.name)
#'
#'
#' save.image(file=paste0(out.dir.name,"re_save_peng.RData"))
#' savehistory(file=paste0(out.dir.name,"re_save_peng.Rhistory"))
#'
ParserBedFile4PengVenn<-function(dir.name,input.file.pattern,out.dir.name){

  file.name=paste0(dir.name,dir(dir.name,recursive = TRUE,pattern=input.file.pattern))
  file.name.2<-as.list(file.name)

  names(file.name.2)=sapply(strsplit(file.name,split="\\/"),"[[",7)

  print(file.name.2)

  #file.name.3<-file.name.2[-6]

  re.out<-lapply(file.name.2,function(u){
    re=read.table(u,header=FALSE)
    re<-as.character(re[,1])
    #colnames(re)=c("Count","GeneName")
    re
  })

  sample.name<-sapply(strsplit(names(re.out),split="_peaks_"),"[[",1)

  names(re.out)=sample.name

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


  venn.plot <- venn.diagram(
    x = re.out[c(1,2)],
    filename = paste0(out.dir.name,names(re.out)[1],names(re.out)[2],"_overlap_venn.tiff"),
    height = 3000,
    width = 3500,
    resolution = 1000,
    col = "black",
    lty = "dotted",
    lwd = 1,
    fill = c("red","blue"),
    alpha = 0.50,
    label.col = c(rep("white",3)),
    cex = 0.5,
    fontfamily = "serif",
    fontface = "bold",
    cat.col = c("red","blue"),
    cat.cex = 0.5,
    cat.pos = 0.5,
    cat.dist = 0.05,
    cat.fontfamily = "serif"
  )

  venn.plot <- venn.diagram(
    x = re.out[c(1,3,5)],
    filename =paste0(out.dir.name,names(re.out)[1],names(re.out)[3],names(re.out)[5],"_overlap_venn.tiff"),
    col = "black",
    lty = "dotted",
    lwd = 2,
    fill = c("red", "orange", "blue"),
    alpha = 0.50,
    label.col = c(rep("white",7)),
    cex = 1,
    fontfamily = "serif",
    fontface = "bold",
    cat.col = c("red", "orange", "blue"),
    cat.cex = 0.8,
    cat.fontfamily = "serif"
  )


  return(re.out)
}
