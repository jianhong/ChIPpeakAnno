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
  
  plotMeta1(sm, xcoords = c(-50000, 50000),overlay=TRUE,profile.names=c(type.index),main=paste(title,"non-NAD to NAD boundary analysis"),lty=3)
  
  
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

# input.bw.path <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/BoundaryAnalysis/InputFiles"
# par(mfrow=c(1,2))

##title="nonXL"
# index <- c(rep(18,2))
# type <- rep(c("G","N"))
# r18 <- BoundaryAnalysis4SingleSample(NAD.peaks[[1]],input.bw.path,type,index,title)

# index <- c(rep(29,2))
# type <- rep(c("G","N"))
# r29 <- BoundaryAnalysis4SingleSample(NAD.peaks[[1]],input.bw.path,type,index,title)

# index <- c(rep(33,2))
# type <- rep(c("G","N"))
# r33 <- BoundaryAnalysis4SingleSample(NAD.peaks[[1]],input.bw.path,type,index,title)

df.nonXL <- data.frame(bases=r18$bases,Group=r18$Group,Coverage=(r18$Coverage+r29$Coverage+r33$Coverage)/3)

label <- title

g <- ggplot(df.nonXL, aes(x=bases, y=Coverage, group=Group)) +
  geom_line(aes(linetype=Group)) + theme(
  legend.position = c(.45, .95),
  legend.justification = c("left", "top"),
  legend.box.just = "left",
  legend.margin = margin(6, 6, 6, 6)
)

g <- g + ggtitle(label) +theme(plot.title=element_text(hjust=0.5),legend.title = element_text(hjust = 0.5,face="bold",size=10))

g <- g + grids(linetype = "dashed")

g

##title="XL"
# index <- c(rep(24,2))
# type <- rep(c("G","N"))
# r24 <- BoundaryAnalysis4SingleSample(NAD.peaks[[2]],input.bw.path,type,index,title)

# index <- c(rep(26,2))
# type <- rep(c("G","N"))
# r26 <- BoundaryAnalysis4SingleSample(NAD.peaks[[2]],input.bw.path,type,index,title)

# index <- c(rep(28,2))
# type <- rep(c("G","N"))
# r28 <- BoundaryAnalysis4SingleSample(NAD.peaks[[2]],input.bw.path,type,index,title)

df.XL <- data.frame(bases=r24$bases,Group=r24$Group,Coverage=(r24$Coverage+r26$Coverage+r28$Coverage)/3)

label <- title

g <- ggplot(df.XL, aes(x=bases, y=Coverage, group=Group)) +
  geom_line(aes(linetype=Group)) + theme(
    legend.position = c(.35, .99),
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.margin = margin(6, 6, 6, 6)
  )

g <- g + ggtitle(label) +theme(plot.title=element_text(hjust=0.5),legend.title = element_text(hjust = 0.5,face="bold",size=10))

g <- g + grids(linetype = "dashed")

g

BoundaryAnalysis4SingleSample <- function(NADpeaks,input.bw.path,type,index,title) {
  
  index.1 <- unique(index)
  type.index <- paste(type,index.1,sep="")
  
  nonNAD.2.NAD <- NADpeaks
  width(nonNAD.2.NAD) <- 1
  nonNAD.2.NAD.region <- nonNAD.2.NAD 
  dd <- 50000
  start(nonNAD.2.NAD.region) <- start(nonNAD.2.NAD) - dd
  end(nonNAD.2.NAD.region) <- end(nonNAD.2.NAD) + dd
  
  files <- dir(input.bw.path, "bw")
  index.nonXL <- c(grep(index.1[1],files))
  targets <- file.path(input.bw.path,files)[index.nonXL]
  
  sm = ScoreMatrixList(targets = targets, windows = nonNAD.2.NAD.region,bin.num = 200)

  re <- plotMeta1(sm, xcoords = c(-50000, 50000),overlay=TRUE,profile.names=c(type.index),main=paste(title,"non-NAD to NAD boundary analysis"))
  
  x  <- re$xcoords
  y1 <- re$metas[[1]]
  y2 <- re$metas[[2]]
  
  NAD.2.nonNAD <- NADpeaks
  
  start(NAD.2.nonNAD) <- end(NAD.2.nonNAD)
  NAD.2.nonNAD.region <- NAD.2.nonNAD
  
  dd <- 50000
  start(NAD.2.nonNAD.region) <- start(NAD.2.nonNAD) - dd
  end(NAD.2.nonNAD.region) <- end(NAD.2.nonNAD) + dd
  
  sm.2 = ScoreMatrixList(targets = targets, windows = NAD.2.nonNAD.region, bin.num = 200)
  
  re.2 <- plotMeta1(sm.2, xcoords = c(-50000, 50000),overlay=TRUE,profile.names=c(type.index),main=paste(title, "NAD to non-NAD boundary analysis",sep=" "))
  
  x  <- re.2$xcoords
  y1r <- rev(re.2$metas[[1]])
  y2r <- rev(re.2$metas[[2]])
  
  y11 <- (y1 + y1r)/2
  y22 <- (y2 + y2r)/2
  
  df <- as.data.frame(rbind(cbind(x,rep("G",length(y11)),y11),cbind(x,rep("N",length(y22)),y22)),stringsAsFactors=F)
  colnames(df) <- c("bases","Group","Coverage")
  df$bases <- as.numeric(df$bases)
  df$Coverage <- as.numeric(df$Coverage)
  
  label <- paste0(title,":",index.1)
  
  g <- ggplot(df, aes(x=bases, y=Coverage, group=Group)) +
    geom_line(aes(linetype=Group))+theme(legend.position="bottom")
  g <- g + ggtitle(label) +theme(plot.title=element_text(hjust=0.5),legend.title = element_text(hjust = 0.5,face="bold",size=10))
  g
  
  df
  
}

output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/BoundaryAnalysis/CoverageBased"

save.image(file = file.path(output.file.dir,"BoundaryAnalysis.RData"))
savehistory(file = file.path(output.file.dir,"BoundaryAnalysis.Rhistory"))
