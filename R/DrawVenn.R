#' Title
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
DrawVenn <- function(devtools, install_github, Vennerable, venneuler, VennDiagram) {
  devtools::install_github("js229/Vennerable");
  library(Vennerable)

  require(venneuler)
  v <- venneuler(c(A=450, B=1800, "A&B"=230))
  plot(v)
  require(VennDiagram)

  AA=1:2547

  a=6025+128
  b=3803+40
  c=19583

  a
  b
  c

  BB=c(1:6153,6154:25736)

  CC=c(6154:25736,25737:29579)

  AA=c(6026:6153,6154:8462,25737:25776,29650:29719)

  Vstem <- Venn(list(ASXL1=AA,SMC1A=BB,RAD21=CC))

  pdf("vennerable.pdf")
  plot(Vstem, doWeights = TRUE,show=list(FaceText="weight",Faces=TRUE))
  #grid.text("Title", vp = viewport(x=0.5, y=.9, w=unit(1, "npc"), h=unit(1, "npc")))

  SetLabels <- VennGetSetLabels(Vstem)

  Cstem3 <- compute.Venn(Vstem,doWeights=TRUE)
  plot(Cstem3)

  SetLabels <- VennGetSetLabels(Cstem3)
  FaceLabels <- VennGetFaceLabels(Cstem3)

  #VennGetSetLabels(Cstem3)

  FaceLabels[FaceLabels$FaceName=="010","x"] <- 10
  FaceLabels[FaceLabels$FaceName=="010","y"] <- 3

  SetLabels[SetLabels$Label=="SMC1A","hjust"] <- "right"
  SetLabels[SetLabels$Label=="SMC1A","x"] <- 110
  SetLabels[SetLabels$Label=="SMC1A","y"] <- 58.78653

  SetLabels[SetLabels$Label=="RAD21","x"] <- -80
  SetLabels[SetLabels$Label=="RAD21","y"] <- 58.78653

  Cstem3 <- VennSetSetLabels(Cstem3,SetLabels)

  Cstem3<-VennSetFaceLabels(Cstem3,FaceLabels)

  grid.newpage()
  plot(Cstem3)


  dev.off()

  dev.off()


  venn.diagram(list(SMC1A=BB,RAD21=CC,ASXL1=AA),fill = c("red", "green","blue"),
               alpha = c(0.5,0.5,0.5), cex = 2,cat.fontface = 4,lty =2,direct.area=TRUE,area.vector=c(25736,23426,2547),fontfamily =3,filename = "trial3.tiff");

  length(intersect(BB,CC))


  Vcombo <- Venn(SetNames = c("SMC1A", "ASXL1", "RAD21"),Weight = c(5387,6025,70,128,3803,17274,40,2309))
  plot(Vcombo)

  plot(Vcombo,type = "squares", show = list(FaceText = "weight",SetLabels = FALSE))
}


