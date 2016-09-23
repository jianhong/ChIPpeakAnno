#' DrawVenn
#'
#' @param devtools
#' @param install_github
#' @param Vennerable
#' @param venneuler
#' @param VennDiagram
#'
#' @return
#' @export
#'
#' @examples
#'
#' out.dir.name="/media/H_driver/2016/Yang/Results/"
#' DrawVenn(out.dir.name)
#'
DrawVenn <- function(out.dir.name) {

  devtools::install_github("js229/Vennerable");
  library(Vennerable)

  BB=c(1:6153,6154:25736)

  CC=c(6154:25736,25737:29579)

  AA=c(6026:6153,6154:8462,25737:25776,29650:29719)

  Vstem <- Venn(list(ASXL1=AA,SMC1A=BB,RAD21=CC))

  SetLabels <- VennGetSetLabels(Vstem)

  Cstem3 <- compute.Venn(Vstem,doWeights=TRUE)

  SetLabels <- VennGetSetLabels(Cstem3)
  FaceLabels <- VennGetFaceLabels(Cstem3)

  FaceLabels[FaceLabels$FaceName=="010","x"] <- 10
  FaceLabels[FaceLabels$FaceName=="010","y"] <- 3

  SetLabels[SetLabels$Label=="SMC1A","hjust"] <- "right"
  SetLabels[SetLabels$Label=="SMC1A","x"] <- 110
  SetLabels[SetLabels$Label=="SMC1A","y"] <- 58.78653

  SetLabels[SetLabels$Label=="RAD21","x"] <- -80
  SetLabels[SetLabels$Label=="RAD21","y"] <- 58.78653

  SetLabels[SetLabels$Label=="ASXL1","y"] <-80

  Cstem3 <- VennSetSetLabels(Cstem3,SetLabels)

  Cstem3<-VennSetFaceLabels(Cstem3,FaceLabels)

  Cstem3@FaceLabels[which(Cstem3@FaceLabels$FaceName=="010"),]$x<-90
  Cstem3@FaceLabels[which(Cstem3@FaceLabels$FaceName=="010"),]$y<--35

  Cstem3@FaceLabels[which(Cstem3@FaceLabels$FaceName=="001"),]$x<--90
  Cstem3@FaceLabels[which(Cstem3@FaceLabels$FaceName=="001"),]$y<--35

  Cstem3@FaceLabels[which(Cstem3@FaceLabels$FaceName=="101"),]$x<--18
  #Cstem3@FaceLabels[which(Cstem3@FaceLabels$FaceName=="101"),]$y<--
  Cstem3@FaceLabels[which(Cstem3@FaceLabels$FaceName=="101"),]$vjust<-"center"

  Cstem3@FaceLabels[which(Cstem3@FaceLabels$FaceName=="110"),]$x<-19
  Cstem3@FaceLabels[which(Cstem3@FaceLabels$FaceName=="110"),]$y<-49.86767

  c.3.set=sum(Cstem3@IndicatorWeight[,4])
  Cstem3@IndicatorWeight[which(row.names(Cstem3@IndicatorWeight)=="000"),4]=35000-c.3.set

  Cstem3@FaceLabels[which(Cstem3@FaceLabels$FaceName=="DarkMatter"),]$y<--130

  temp<-ol$venn_cnt

  index.A=grep(colnames(Cstem3@IndicatorWeight)[1],colnames(temp))
  index.B=grep("SMC1",colnames(temp))
  index.C=grep("Rad21",colnames(temp))

  temp2<-as.data.frame(temp[,c(index.A,index.B,index.C,4)])
  colnames(temp2)[1:3]=colnames(Cstem3@IndicatorWeight)[1:3]

  row.names(temp2)<-with(temp2,paste0(temp2$ASXL1,temp2$SMC1A,temp2$RAD21))

  row.names(temp2)<-row.names(Cstem3@IndicatorWeight)

  temp3<-merge(Cstem3@IndicatorWeight,temp2,by=0,sort=FALSE)

  Cstem3@IndicatorWeight[,4]<-temp3$Counts

  c.3.set=sum(Cstem3@IndicatorWeight[,4])
  Cstem3@IndicatorWeight[which(row.names(Cstem3@IndicatorWeight)=="000"),4]=35000-c.3.set

  grid.newpage()

  pdf(paste0(out.dir.name,"venn.pdf"))
  plot(Cstem3)
  dev.off()

}


