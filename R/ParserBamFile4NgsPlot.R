#' ParserBamFile4NgsPlot
#'
#' @param dir.name
#' @param input.file.pattern
#'
#' @return
#' @export
#'
#' @examples
#'

#'
#' sh ~/Code/BashRunMACS1-4-2_4_Danny_chip_seq3.sh /scratch/projects/bbc/aiminy_project/Bam_sorted/ 7 ".sorted.bam$" /scratch/projects/bbc/aiminy_project/Bam_marked_sorted/
#'
#' dir.name="/scratch/projects/bbc/aiminy_project/Bam_sorted/"
#' index=7
#' input.file.pattern=".sorted.bam$"
#' out.dir.name="/scratch/projects/bbc/aiminy_project/Bam_marked_sorted/"
#' outfile="test.txt"
#'
#' re.from.bed.peng.4.venn<-ParserBamFile4NgsPlot(dir.name,index,input.file.pattern,out.dir.name,outfile)
#'
#' save.image(file=paste0(out.dir.name,"re_save_peng.RData"))
#' savehistory(file=paste0(out.dir.name,"re_save_peng.Rhistory"))
#'
ParserBamFile4NgsPlot<-function(dir.name,index,input.file.pattern,out.dir.name,outfile){

  file.name=paste0(dir.name,dir(dir.name,recursive = TRUE,pattern=input.file.pattern))
  file.name.2<-as.list(file.name)


  #file.name="/scratch/projects/bbc/Project/Danny_chip/Alignment/BWA/2016-10-14-Hyunho-Yoon-02-MDA-MB-231-2-cJun_S5_R1.marked.bam"

  #strsplit(file.name,split="\\/")

  names(file.name.2)=sapply(strsplit(file.name,split="\\/"),"[[",index)

  print(file.name.2)

  ab<-names(file.name.2)

  #ab<-"2016-10-14-Hyunho-Yoon-13-MDA-MB-231-DD-2-p27_S15_R1.marked.bam"

  ab2<-sapply(strsplit(ab,split="\\."),"[[",1)

  print(ab2)
  #ab3<-sapply(strsplit(ab2,split="\\-"),"[[",length(strsplit(ab2,split="\\-")[[1]]))


  re<-cbind(paste0(dir.name,names(file.name.2)),rep(-1,length(names(file.name.2))),ab2)

  write.table(re, file =paste0(out.dir.name,outfile), append = FALSE, quote = F, sep = " ",eol = "\n", na = "NA", dec = ".", row.names = F, col.names = F)


  cmd1="ngs.plot.r -G hg19 -R tss -C"
  cmd2="-O"
  cmd3="-L 4000"
  #cmd3="-L 4000 -RR 1 -CD 1 -CO \\\"blue\\\""

  #ngs.plot.r -G hg19 -R tss -C $1 -O $2 -L 4000 -RR 1 -CD 1 -CO "blue"


  #file.name.3<-file.name.2[-6]

  re.out<-lapply(file.name.2,function(u){

     cmd4=paste(cmd1,u,cmd2,u,cmd3,sep=" ")

     print(cmd4)
     system(cmd4, intern = TRUE, ignore.stderr = TRUE)

     #re=read.table(u,header=FALSE)
  #  re<-as.character(re[,1])
  #  #colnames(re)=c("Count","GeneName")
  #  re
  })

  #sample.name<-sapply(strsplit(names(re.out),split="_peaks_"),"[[",1)

  #names(re.out)=sample.name

  # A=re.out[[1]] #10_2470IUPUI_WT_BM_ASXL1
  # B=re.out[[2]] #10_WT_BM_ASXL1
  # C=re.out[[3]] #11_2470IUPUI_WT_BM_SMC1
  # D=re.out[[5]] #13_2470IUPUI_WT_BM_Rad21
  #
  # A.name=sample.name[1] #10_2470IUPUI_WT_BM_ASXL1
  # B.name=sample.name[2] #10_WT_BM_ASXL1
  # C.name=sample.name[3] #11_2470IUPUI_WT_BM_SMC1
  # D.name=sample.name[5] #13_2470IUPUI_WT_BM_Rad21
  #
  # #A <- createRandomRegions(nregions=50, length.mean=5000000, length.sd=3000000)
  # #B <- c(A[1:25], createRandomRegions(nregions=25, length.mean=500000, length.sd=30000))
  # #class(A)
  # #class(B)
  #
  #
  # # TestOverlap <- function(A, B, A.name,B.name,out.dir.name,n_permutations) {
  # #   NumOverLap<-numOverlaps(A, B, count.once=TRUE)
  # #   pt <- overlapPermTest(A=A, B=B, ntimes=n_permutations)
  # #   png(filename=paste0(out.dir.name,A.name,B.name,n_permutations,"_overlap_test.png"))
  # #   plot(pt)
  # #   dev.off()
  # #
  # #   lz <- localZScore(pt=pt, A=A, B=B)
  # #   png(filename=paste0(out.dir.name,A.name,B.name,n_permutations,"_localZScore.png"))
  # #   plot(lz)
  # #   dev.off()
  # # }
  # #
  # # TestOverlap(A,B,A.name,B.name,out.dir.name,100)
  # # TestOverlap(A,C,A.name,C.name,out.dir.name,100)
  # # TestOverlap(A,D,A.name,D.name,out.dir.name,100)
  # #
  # # TestOverlap(B,C,B.name,C.name,out.dir.name,100)
  # # TestOverlap(B,D,B.name,D.name,out.dir.name,100)
  # # TestOverlap(C,D,C.name,D.name,out.dir.name,100)
  # #
  #
  # # temp.name<-strsplit(names(file.name.2),split="\\.")
  #
  # #a = readBed(system.file("extdata", "examples/combined_regions.bed",
  # #                        package="LOLA"))
  #
  #
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


  #return(re.out)
}
