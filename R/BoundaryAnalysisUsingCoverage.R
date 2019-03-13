# input.bw.path <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/BoundaryAnalysis/InputFiles"
# par(mfrow=c(3,2))
# 
# index <- c(rep(18,2),rep(29,2),rep(33,2))
# type <- rep(c("G","N"),3)
# title="nonXL"
# BoundaryAnalysis2(NAD.peaks[[1]],input.bw.path,type,index,title)
# 
# index <- c(rep(24,2),rep(26,2),rep(28,2))
# type <- rep(c("G","N"),3)
# title="XL"
# BoundaryAnalysis2(NAD.peaks[[2]],input.bw.path,type,index,title)

BoundaryAnalysis2 <- function(NADpeaks,input.bw.path,type,index,title) {
  
  index.1 <- unique(index)
  type.index <- paste(type,index,sep="")
  
  #NADpeaks <- NAD.peaks[[1]]
  
  nonNAD.2.NAD <- NADpeaks
  width(nonNAD.2.NAD) <- 1
  nonNAD.2.NAD.region <- nonNAD.2.NAD 
  dd <- 50000
  start(nonNAD.2.NAD.region) <- start(nonNAD.2.NAD) - dd
  end(nonNAD.2.NAD.region) <- end(nonNAD.2.NAD) + dd
  
  files <- dir(input.bw.path, "bw")
  index.nonXL <- c(grep(index.1[1],files),grep(index.1[2],files),grep(index.1[3],files))
  targets <- file.path(input.bw.path,files)[index.nonXL]
  
  sm = ScoreMatrixList(targets = targets, windows = nonNAD.2.NAD.region, bin.num = 200)
  
  title <- "nonXL"
  
  plotMeta(sm, xcoords = c(-50000, 50000),overlay=TRUE,profile.names=c(type.index),main=paste(title,"non-NAD to NAD boundary analysis"))
  
  
  NAD.2.nonNAD <- NADpeaks
  
  start(NAD.2.nonNAD) <- end(NAD.2.nonNAD)
  
  NAD.2.nonNAD.region <- NAD.2.nonNAD
  
  dd <- 50000
  start(NAD.2.nonNAD.region) <- start(NAD.2.nonNAD) - dd
  end(NAD.2.nonNAD.region) <- end(NAD.2.nonNAD) + dd
  
  #targets1 <- "~/Aimin/DropboxUmass/NADfinder/Aimin/BoundaryAnalysis/N18.sort.markDup.bw"
  
  sm = ScoreMatrixList(targets = targets, windows = NAD.2.nonNAD.region, bin.num = 200)
  
  plotMeta(sm, xcoords = c(-50000, 50000),overlay=TRUE,profile.names=c(type.index),main=paste(title, "NAD to non-NAD boundary analysis",sep=" "))
  
}
