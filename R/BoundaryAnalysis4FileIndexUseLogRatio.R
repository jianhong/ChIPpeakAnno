#!/usr/bin/env Rscript

BoundaryAnalysis4FileIndexUseLogRatio <- function(NAD.peaks.F121.9, input.bw.path, title, file.index, rep.index){
  
  nonNAD.2.NAD <- NAD.peaks.F121.9[[1]]
  width(nonNAD.2.NAD) <- 1
  nonNAD.2.NAD.region <- nonNAD.2.NAD 
  dd <- 50000
  start(nonNAD.2.NAD.region) <- start(nonNAD.2.NAD) - dd
  end(nonNAD.2.NAD.region) <- end(nonNAD.2.NAD) + dd
  
  files <- dir(input.bw.path)
  index.nonXL <- c(file.index)
  targets <- file.path(input.bw.path,files)[index.nonXL]
  
  sm.1 = ScoreMatrixBin(target = targets[1], windows = nonNAD.2.NAD.region,bin.num = 200)
  sm.2 = ScoreMatrixBin(target = targets[2], windows = nonNAD.2.NAD.region,bin.num = 200)
  
  sm.2.1 <- (sm.2+1)/(sm.1+1)
  
  sm.log <- log2(sm.2.1)
  
  re <- plotMeta1(sm.log,xcoords = c(-50000, 50000),overlay=FALSE,main=paste(title,"non-NAD to NAD boundary analysis"))
  
  x  <- re$xcoords
  y1 <- re$metas[[1]]
  
  NAD.2.nonNAD <- NAD.peaks.F121.9[[1]]
  start(NAD.2.nonNAD) <- end(NAD.2.nonNAD)
  NAD.2.nonNAD.region <- NAD.2.nonNAD
  
  dd <- 50000
  start(NAD.2.nonNAD.region) <- start(NAD.2.nonNAD) - dd
  end(NAD.2.nonNAD.region) <- end(NAD.2.nonNAD) + dd
  
  sm.1 = ScoreMatrixBin(target = targets[1], windows = NAD.2.nonNAD.region,bin.num = 200)
  sm.2 = ScoreMatrixBin(target = targets[2], windows = NAD.2.nonNAD.region,bin.num = 200)
  
  sm.2.1 <- (sm.2+1)/(sm.1+1)
  
  sm.log <- log2(sm.2.1)
  
  re.2 <- plotMeta1(sm.log,xcoords = c(-50000, 50000),overlay=FALSE,main=paste(title,"NAD to non-NAD boundary analysis"))
  
  x  <- re.2$xcoords
  y1r <- rev(re.2$metas[[1]])
  
  y11 <- (y1 + y1r)/2
  
  df <- as.data.frame(cbind(x,rep("log2(N/G)",length(y11)),y11),stringsAsFactors=F)
  
  colnames(df) <- c("bases","Group","Coverage")
  df$bases <- as.numeric(df$bases)
  df$Coverage <- as.numeric(df$Coverage)
  
  label <- paste0(title,":",rep.index)
  
  g <- ggplot(df, aes(x=bases, y=Coverage, group=Group)) +
    geom_line(aes(linetype=Group))+theme(legend.position="bottom")
  g <- g + ggtitle(label) +theme(plot.title=element_text(hjust=0.5),legend.title = element_text(hjust = 0.5,face="bold",size=10))
  g
  
  df
  
}

file.index <- c(1,13)
rep.index <- 1
re.1 <- BoundaryAnalysis4FileIndexUseLogRatio(NAD.peaks.F121.9, input.bw.path, title, file.index, rep.index)

file.index <- c(4,16)
rep.index <- 2
re.2 <- BoundaryAnalysis4FileIndexUseLogRatio(NAD.peaks.F121.9,input.bw.path,title, file.index, rep.index)

df.nonXL <- data.frame(bases=re.1$bases,Group=re.1$Group,Coverage=(re.1$Coverage+re.2$Coverage)/2)

title <- "F121.9"
label <- title

g <- ggplot(df.nonXL, aes(x=bases, y=Coverage, group=Group)) +
  geom_line(aes(linetype=Group)) + theme(
    legend.position = c(.45, .99),
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.margin = margin(6, 6, 6, 6)
  ) + ylab("log2Ratio")

g <- g + ggtitle(label) +theme(plot.title=element_text(hjust=0.5),legend.title = element_text(hjust = 0.5,face="bold",size=10))

g <- g + grids(linetype = "dashed")

g <- g + geom_hline(yintercept=0, linetype="dashed", color = "red")

g

ggsave(file.path(output.file.dir,"F121_9_logRatio.png"))

save.image(file = file.path(output.file.dir,"Boundary4_logRatio_N_G.RData"))
savehistory(file = file.path(output.file.dir,"Boundary4_logRatio_N_G.Rhistory"))
