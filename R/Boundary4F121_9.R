#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("genomation")
BiocManager::install("ComplexHeatmap")
install.packages("ggpubr")

library(genomation)
library(ComplexHeatmap)
library(ChIPpeakAnno)
library(ggplot2)
library(grid)
library(ggpubr)

F121.9.NAD <- '/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/avesig filtered NADs/F121_9_all_1.5FC_F200_minp_adjpVal_avesig_greater_than_1.7_NAD_greater_50kb.bed'
x <- F121.9.NAD
peaks=read.table(x,header=F,skip=1)
colnames(peaks)[1:4]= c("chr","start","end","pValue")
peaks=toGRanges(peaks)
NAD.peaks.F121.9 <- list(F121.9.NAD=peaks)

output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/BoundaryAnalysis4F121_9"

input.bw.path <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/Aizhan.allele.specific/DNASEQ/hybrid_align/markDup_align/batch1_bw"

# Notes
#replicate 1 nucleolar DNA:
#  Dropbox/pub_old_galaxy/F121-9/Aizhan allele-specific/DNASEQ/hybrid_align/markDup_align/batch1_bw/nusDNA_1.markDup.bw
#replicate 1 genomic DNA:
#  Dropbox/pub_old_galaxy/F121-9/Aizhan allele-specific/DNASEQ/hybrid_align/markDup_align/batch1_bw/gDNA_1.markDup.bw

#replicate 2 nucleolar DNA:
#  Dropbox/pub_old_galaxy/F121-9/Aizhan allele-specific/DNASEQ/hybrid_align/markDup_align/batch1_bw/nusDNA_2.markDup.bw
#replicate 2 genomic DNA:
#  Dropbox/pub_old_galaxy/F121-9/Aizhan allele-specific/DNASEQ/hybrid_align/markDup_align/batch1_bw/gDNA_2.markDup.bw

BoundaryAnalysis4FileIndex <- function(NAD.peaks.F121.9, input.bw.path, title, file.index, rep.index) {
  
  nonNAD.2.NAD <- NAD.peaks.F121.9[[1]]
  width(nonNAD.2.NAD) <- 1
  nonNAD.2.NAD.region <- nonNAD.2.NAD 
  dd <- 50000
  start(nonNAD.2.NAD.region) <- start(nonNAD.2.NAD) - dd
  end(nonNAD.2.NAD.region) <- end(nonNAD.2.NAD) + dd
  
  #files <- dir(input.bw.path, "bw")
  files <- dir(input.bw.path)
  index.nonXL <- c(file.index)
  targets <- file.path(input.bw.path,files)[index.nonXL]
  
  sm = ScoreMatrixList(targets = targets, windows = nonNAD.2.NAD.region,bin.num = 200)
  
  re <- plotMeta1(sm, xcoords = c(-50000, 50000),overlay=TRUE,profile.names=c("G1","N1"),main=paste(title,"non-NAD to NAD boundary analysis"))
  
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
  
  re.2 <- plotMeta1(sm.2, xcoords = c(-50000, 50000),overlay=TRUE,profile.names=c("G1","N1"),main=paste(title, "NAD to non-NAD boundary analysis",sep=" "))
  
  
  x  <- re.2$xcoords
  y1r <- rev(re.2$metas[[1]])
  y2r <- rev(re.2$metas[[2]])
  
  y11 <- (y1 + y1r)/2
  y22 <- (y2 + y2r)/2
  
  df <- as.data.frame(rbind(cbind(x,rep("G",length(y11)),y11),cbind(x,rep("N",length(y22)),y22)),stringsAsFactors=F)
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
re.1 <- BoundaryAnalysis4FileIndex(NAD.peaks.F121.9, input.bw.path, title, file.index, rep.index)

file.index <- c(4,16)
rep.index <- 2
re.2 <- BoundaryAnalysis4FileIndex(NAD.peaks.F121.9,input.bw.path,title, file.index, rep.index)

df.nonXL <- data.frame(bases=re.1$bases,Group=re.1$Group,Coverage=(re.1$Coverage+re.2$Coverage)/2)

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

input.bw.path.mef <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/BoundaryAnalysis/InputFiles"
files.mef <- dir(input.bw.path.mef)

file.index <- c(1,8)
rep.index <- 1
re.11 <- BoundaryAnalysis4FileIndex(NAD.peaks.F121.9,input.bw.path.mef,title, file.index, rep.index)

file.index <- c(6,13)
rep.index <- 2
re.22 <- BoundaryAnalysis4FileIndex(NAD.peaks.F121.9,input.bw.path.mef,title, file.index, rep.index)

file.index <- c(7,14)
rep.index <- 3
re.33 <- BoundaryAnalysis4FileIndex(NAD.peaks.F121.9,input.bw.path.mef,title, file.index, rep.index)


df.nonXL <- data.frame(bases=re.11$bases,Group=re.11$Group,Coverage=(re.11$Coverage+re.22$Coverage+re.33$Coverage)/3)

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

NAD.peaks[[1]]
files <- dir(input.bw.path)

file.index <- c(1,8)
rep.index <- 1
re.11 <- BoundaryAnalysis4FileIndex(NAD.peaks,input.bw.path,title, file.index, rep.index)

file.index <- c(6,13)
rep.index <- 2
re.22 <- BoundaryAnalysis4FileIndex(NAD.peaks,input.bw.path,title, file.index, rep.index)

file.index <- c(7,14)
rep.index <- 3
re.33 <- BoundaryAnalysis4FileIndex(NAD.peaks,input.bw.path,title, file.index, rep.index)

df.nonXL <- data.frame(bases=re.11$bases,Group=re.11$Group,Coverage=(re.11$Coverage+re.22$Coverage+re.33$Coverage)/3)

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


unfiltered.F121.NAD <- '/Users/aiminyan/Aimin/DropboxUmass/NADfinder/BedFiles/Aizhan F121-9 comparisons/F121_9_all_1.5FC_F200_minp_adjPval copy.bed'

peaks=read.table(unfiltered.F121.NAD,header=F,skip=1)
colnames(peaks)[1:4]= c("chr","start","end","pValue")
peaks=toGRanges(peaks)
NAD.peaks.F121.9 <- list(F121.9.NAD=peaks)


input.bw.path <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/Aizhan.allele.specific/DNASEQ/hybrid_align/markDup_align/batch1_bw"

file.index <- c(1,13)
rep.index <- 1
re.1 <- BoundaryAnalysis4FileIndex(NAD.peaks.F121.9, input.bw.path, title, file.index, rep.index)

file.index <- c(4,16)
rep.index <- 2
re.2 <- BoundaryAnalysis4FileIndex(NAD.peaks.F121.9,input.bw.path,title, file.index, rep.index)

df.nonXL <- data.frame(bases=re.1$bases,Group=re.1$Group,Coverage=(re.1$Coverage+re.2$Coverage)/2)

title <- "unfiltered.F121.9"
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

save.image(file = file.path(output.file.dir,"Boundary4_unfilteredF121_9.RData"))
savehistory(file = file.path(output.file.dir,"Boundary4_unfilteredF121_9.Rhistory"))

save.image(file = file.path(output.file.dir,"BoundaryAnalysis4_F121_9.RData"))
savehistory(file = file.path(output.file.dir,"BoundaryAnalysis4_F121_9.Rhistory"))

save.image(file = file.path(output.file.dir,"BoundaryAnalysis4_F121_9_new_function.RData"))
savehistory(file = file.path(output.file.dir,"BoundaryAnalysis4_F121_9_new_function.Rhistory"))

save.image(file = file.path(output.file.dir,"BoundaryAnalysis4_F121_9_Fig_E.RData"))
savehistory(file = file.path(output.file.dir,"BoundaryAnalysis4_F121_9_Fig_E.Rhistory"))

