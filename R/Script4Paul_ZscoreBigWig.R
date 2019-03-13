# sample.name <- "N18"
# getBigWig2(sample.name,range.rato.ss,output.dir)

# sample.name <- "N29"
# getBigWig2(sample.name,range.rato.ss,output.dir)

# sample.name <- "N33"
# getBigWig2(sample.name,range.rato.ss,output.dir)  

getBigWig2 <- function(sample.name,range.rato.ss,output.dir) {
  
  range.rato.ss$names <- tools::file_path_sans_ext(range.rato.ss$names)
  
  s.index <- grep(sample.name,colnames(range.rato.ss))
  df.2.gr.temp <- range.rato.ss[,c(4,1,2,s.index)]
  colnames(df.2.gr.temp) <- c("seqnames","start","end","score")
  
  df.2.gr <-  GRanges(df.2.gr.temp)
  seqlengths(df.2.gr) <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)[match(names(seqlengths(df.2.gr)),names(seqlengths(BSgenome.Mmusculus.UCSC.mm10)))] 
  
  len <- unlist(lapply(seList, length))
  ZZ <- lapply(1:length(len), function(u,len,df.2.gr){
    x <- names(len)[u]
    y <- len[u]
    z <- df.2.gr[which(seqnames(df.2.gr)==x)]
    non.olp.index <- getNonOverLappingPeakIndex(z)
    z <- z[non.olp.index]
    z
  },len,df.2.gr)
  ZZZ <- do.call(base::c,ZZ)
  
  if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}
  test_bw_out <- file.path(output.dir, paste0(tools::file_path_sans_ext(colnames(range.rato.ss)[s.index]),".bw"))
  export.bw(ZZZ,test_bw_out)
}

# df.2.gr.temp <- range.rato.ss[,c(4,1,2,5)]
# colnames(df.2.gr.temp) <- c("seqnames","start","end","score")
# 
# df.2.gr <-  GRanges(df.2.gr.temp)
# seqlengths(df.2.gr) <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)[match(names(seqlengths(df.2.gr)),names(seqlengths(BSgenome.Mmusculus.UCSC.mm10)))] 
# 
# output.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/BoundaryAnalysis"
# 
# if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}
# test_bw_out <- file.path(output.dir, paste0(tools::file_path_sans_ext(colnames(range.rato.ss)[s.index]),".bw"))
# export.bw(ZZZ,test_bw_out)

par(mfrow=c(2,1))
input.bw.path <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/BoundaryAnalysis"
index <- c(rep(18,1),rep(29,1),rep(33,1))
type <- rep(c("N"),1)
title="nonXL"
BoundaryAnalysis3(NAD.peaks[[1]],input.bw.path,type,index,title)

BoundaryAnalysis3 <- function(NADpeaks,input.bw.path,type,index,title) {
  
  index.1 <- unique(index)
  type.index <- paste(type,index,sep="")
  
  #NADpeaks <- NAD.peaks[[1]]
  
  NADpeaks <- NADpeaks
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
  
  #title <- "nonXL"
  
  plotMeta(sm, xcoords = c(-50000, 50000),overlay=TRUE,profile.names=c(type.index),main=paste(title,"non-NAD to NAD"))
  
  NAD.2.nonNAD <- NADpeaks
  
  start(NAD.2.nonNAD) <- end(NAD.2.nonNAD)
  
  NAD.2.nonNAD.region <- NAD.2.nonNAD
  
  dd <- 50000
  start(NAD.2.nonNAD.region) <- start(NAD.2.nonNAD) - dd
  end(NAD.2.nonNAD.region) <- end(NAD.2.nonNAD) + dd
  
  #targets1 <- "~/Aimin/DropboxUmass/NADfinder/Aimin/BoundaryAnalysis/N18.sort.markDup.bw"
  
  sm = ScoreMatrixList(targets = targets, windows = NAD.2.nonNAD.region, bin.num = 200)
  
  plotMeta(sm, xcoords = c(-50000, 50000),overlay=TRUE,profile.names=c(type.index),main=paste(title, "NAD to non-NAD",sep=" "))
  
}

range.rato.s <- lapply(1:length(seList), function(u,seList){
  range.ratio <- cbind(as.data.frame(ranges(seList[[u]])),assays(seList[[u]])$bcRatio)
  range.ratio
},seList)
range.rato.ss <- do.call(rbind,range.rato.s)

# sample.name <- "N24"
# getBigWig2(sample.name,range.rato.ss,output.dir)

# sample.name <- "N26"
# getBigWig2(sample.name,range.rato.ss,output.dir)

# sample.name <- "N28"
# getBigWig2(sample.name,range.rato.ss,output.dir)  

par(mfrow=c(2,1))
input.bw.path <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/BoundaryAnalysis"
index <- c(rep(24,1),rep(26,1),rep(28,1))
type <- rep(c("N"),1)
title="XL"
BoundaryAnalysis3(NAD.peaks[[2]],input.bw.path,type,index,title)

range.rato.z <- lapply(1:length(seList), function(u,seList){
  range.ratio <- cbind(as.data.frame(ranges(seList[[u]])),assays(seList[[u]])$zscore)
  range.ratio
},seList)

range.rato.zz <- do.call(rbind,range.rato.z)

output.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/BoundaryAnalysis/ZscoreBased"

# sample.name <- "N24"
# getBigWig2(sample.name,range.rato.zz,output.dir)

# sample.name <- "N26"
# getBigWig2(sample.name,range.rato.zz,output.dir)

# sample.name <- "N28"
# getBigWig2(sample.name,range.rato.zz,output.dir)  

par(mfrow=c(2,2))
input.bw.path <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/BoundaryAnalysis/ZscoreBased"
index <- c(rep(24,1),rep(26,1),rep(28,1))
type <- rep(c("N"),1)
title="XL"
y_title="zscore"
BoundaryAnalysis3(NAD.peaks[[2]],input.bw.path,type,index,title,y_title=y_title)

#/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output/Aizhan_3_sets_bed_4/XL.only.bed
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

range.rato.z1 <- lapply(1:length(seList), function(u,seList){
  
  range01(apply(assays(seList[[u]])$zscore,1,mean))
  
  range.ratio <- cbind(as.data.frame(ranges(seList[[u]])),assays(seList[[u]])$zscore)
  range.ratio
},seList)

range.rato.zz1 <- do.call(rbind,range.rato.z1)

sample.name <- "N18"
getBigWig2(sample.name,range.rato.zz1,output.dir)

sample.name <- "N29"
getBigWig2(sample.name,range.rato.zz1,output.dir)

sample.name <- "N33"
getBigWig2(sample.name,range.rato.zz1,output.dir)

input.bw.path <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/BoundaryAnalysis/ZscoreBased"
index <- c(rep(18,1),rep(29,1),rep(33,1))
type <- rep(c("N"),1)
title="nonXL"
BoundaryAnalysis3(NAD.peaks[[1]],input.bw.path,type,index,title,y_title=y_title)

range.rato.z2 <- lapply(1:length(seList), function(u,seList){
  nonXL = range01(apply(assays(seList[[u]])$zscore,1,mean))
  range.ratio <- cbind(as.data.frame(ranges(seList[[u]])),nonXL=nonXL)
  range.ratio
},seList)

range.rato.zz2 <- do.call(rbind,range.rato.z2)
                          
sample.name <- "nonXL"
getBigWig2(sample.name,range.rato.zz2,output.dir)

BoundaryAnalysis4 <- function(NADpeaks,input.bw.file,y_title) {
  
  NADpeaks <- NADpeaks
  nonNAD.2.NAD <- NADpeaks
  width(nonNAD.2.NAD) <- 1
  nonNAD.2.NAD.region <- nonNAD.2.NAD 
  dd <- 50000
  start(nonNAD.2.NAD.region) <- start(nonNAD.2.NAD) - dd
  end(nonNAD.2.NAD.region) <- end(nonNAD.2.NAD) + dd
  
  targets <- input.bw.file
  
  sm = ScoreMatrixList(targets = targets, windows = nonNAD.2.NAD.region, bin.num = 200)
  
  plotMeta(sm, xcoords = c(-50000, 50000),overlay=TRUE,profile.names=c(type.index),main=paste(title,"non-NAD to NAD"),ylab=y_title)
  
  NAD.2.nonNAD <- NADpeaks
  
  start(NAD.2.nonNAD) <- end(NAD.2.nonNAD)
  
  NAD.2.nonNAD.region <- NAD.2.nonNAD
  
  dd <- 50000
  start(NAD.2.nonNAD.region) <- start(NAD.2.nonNAD) - dd
  end(NAD.2.nonNAD.region) <- end(NAD.2.nonNAD) + dd
  
  sm = ScoreMatrixList(targets = targets, windows = NAD.2.nonNAD.region, bin.num = 200)
  
  plotMeta(sm, xcoords = c(-50000, 50000),overlay=TRUE,profile.names=c(type.index),main=paste(title, "NAD to non-NAD",sep=" "),ylab=y_title)
  
}

par(mfrow=c(2,2))

input.bw.file <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/BoundaryAnalysis/ZscoreBased/nonXL.bw"
title="nonXL"
y_title="Average zscore"
BoundaryAnalysis4(NAD.peaks[[1]],input.bw.file,y_title=y_title)

range.rato.XL <- lapply(1:length(seList.XL), function(u,seList.XL){
  XL = range01(apply(assays(seList.XL[[u]])$zscore,1,mean))
  range.ratio <- cbind(as.data.frame(ranges(seList.XL[[u]])),XL=XL)
  range.ratio
},seList.XL)

range.rato.XL.L <- do.call(rbind,range.rato.XL)

sample.name <- "XL"
getBigWig2(sample.name,range.rato.XL.L,output.dir)

input.bw.file <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/BoundaryAnalysis/ZscoreBased/XL.bw"
title="XL"
y_title="Average zscore"
BoundaryAnalysis4(NAD.peaks[[2]],input.bw.file,y_title=y_title)

save.image(file = file.path(output.dir,"ZscoreBased.RData"))
savehistory(file = file.path(output.dir,"ZscoreBased.Rhistory"))