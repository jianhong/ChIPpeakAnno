#!/usr/bin/env Rscript

library(ComplexHeatmap)
library(gplots)
library(cowplot)

# input.bed.dir <- "/Users/aiminyan/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/bedFilesforDensityPlots/byClass"

# input.bw.path <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/H3K27me3_bigwig"

# output.file.dir <- "~/Aimin/umw_nathan_lawson/Aimin/H3K27me3_cvg_3"
# cvglists.l <- overLapWithOtherFeatures(input.file.dir,input.bw.path,output.file.dir)
# 
# output.file.dir <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/Aimin/DensityPlot"

# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/test"

# input.bed.dir <- "/Users/aiminyan/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/bedFilesforDensityPlots/byClass"
# input.bw.path <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/H3K27me3_bigwig"
# dd <- 3000
# dd <- 10000
# tfType <- "H3K27me3"
# class_pattern <- "tss"
# tssType <- c("Class1","Class2","Class3","Endo")
# output.file.dir <- paste0("~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/H3K27me3_",class_pattern)
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType)

# dd <- 3000
# dd <- 10000
# class_pattern <- "enh"
# output.file.dir <- paste0("~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/H3K27me3_",class_pattern)
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern)

# dd <- 3000
# dd <- 10000
# class_pattern <- "DHS"
# output.file.dir <- paste0("~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/H3K27me3_",class_pattern)
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern)


# ArteryTSS.bed
# VeinTSS.bed
# CommonTSS.bed
# ArteryEnh.bed
# VeinEnh.bed
# CommonEnh.bed

# input.bed.dir <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/bedFilesforDensityPlots/AllElementsByCellType/"
# input.bw.path <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/H3K27me3_bigwig"
# dd <- 3000
# dd <- 10000
# class_pattern <- "TSS"
# tssType <- c("Artery","Common","Vein")
# output.file.dir <- paste0("~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/H3K27me3_",paste0("Big_",class_pattern))
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType)

# dd <- 3000
# dd <- 10000
# class_pattern <- "Enh"
# tssType <- c("Artery","Common","Vein")
# output.file.dir <- paste0("~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/H3K27me3_",paste0("Big_",class_pattern))
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType)


# input.bed.dir <- "/Users/aiminyan/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/bedFilesforDensityPlots/byClass"
# input.bw.path <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/NR2F2kd_bigwig"
# dd <- 3000
# dd <- 10000
# tfType <- "TRC-NS"
# class_pattern <- "tss"
# tssType <- c("Class1","Class2","Class3","Endo")
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/NR2F2"
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

# dd <- 3000
# dd <- 10000
# tfType <- "TRC-NS"
# class_pattern <- "enh"
# tssType <- c("Class1","Class2","Class3","Endo")
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/NR2F2"
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

# dd <- 3000
# dd <- 10000
# tfType <- "TRC-NS"
# class_pattern <- "DHS"
# tssType <- c("Class1","Class2","Class3","Endo")
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/NR2F2"
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

# dd <- 3000
# dd <- 10000
# tfType <- "shNR2F2"
# class_pattern <- "enh"
# tssType <- c("Class1","Class2","Class3","Endo")
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/NR2F2"
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

# dd <- 3000
# dd <- 10000
# tfType <- "shNR2F2"
# class_pattern <- "DHS"
# tssType <- c("Class1","Class2","Class3","Endo")
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/NR2F2"
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

# tfType <- "shNR2F2"
# class_pattern <- "tss"
# tssType <- c("Class1","Class2","Class3","Endo")
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/NR2F2"
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

# input.bed.dir <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/bedFilesforDensityPlots/AllElementsByCellType/"
# input.bw.path <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/NR2F2kd_bigwig"
# dd <- 3000
# dd <- 10000
# tfType <- "TRC-NS"

# class_pattern <- "TSS"
# tssType <- c("Artery","Common","Vein")
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/NR2F2"
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

# class_pattern <- "Enh"
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/NR2F2"
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

# dd <- 3000
# dd <- 10000
# tfType <- "shNR2F2"

# class_pattern <- "TSS"
# tssType <- c("Artery","Common","Vein")
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/NR2F2"
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

# class_pattern <- "Enh"
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/NR2F2"
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

# input.bed.dir <- "/Users/aiminyan/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/bedFilesforDensityPlots/byClass"
# input.bw.path <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/Aimin/EZH2bam_ENCODE_bw"
# 
# dd <- 3000
# dd <- 10000
# tfType <- "EZH2"
# class_pattern <- "tss"
# class_pattern <- "enh"
# class_pattern <- "DHS"
# tssType <- c("Class1","Class2","Class3","Endo")


# input.bed.dir <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/bedFilesforDensityPlots/AllElementsByCellType/"
#
# class_pattern <- "TSS"
# tssType <- c("Artery","Common","Vein")

#
# class_pattern <- "Enh"
# tssType <- c("Artery","Common","Vein")

# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/EZH2"
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

# input.bed.dir <- "/Users/aiminyan/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/bedFilesforDensityPlots/byClass"
# input.bw.path <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/Aimin/HDAC1cr_bw"
# 
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/HDAC1"
# dd <- 3000
# dd <- 10000
# tfType <- "HDAC1"
# class_pattern <- "tss"
# class_pattern <- "enh"
# class_pattern <- "DHS"
# tssType <- c("Class1","Class2","Class3","Endo")
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/HDAC2"
# input.bed.dir <- "~/Aimin/project/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/bedFilesforDensityPlots/AllElementsByCellType"
# input.bw.path <- "~/Aimin/project/umw_nathan_lawson/Aimin/HDAC2cr_bw"
# class_pattern <- "TSS"
# tssType <- c("Artery","Common","Vein")
# dd <- 3000
# tfType <- "HDAC2"
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)
## dd <- 10000
## null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

# class_pattern <- "Enh"
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)
# dd <- 3000
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)


# input.bed.dir <- "/Users/aiminyan/Aimin/project/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/bedFilesforDensityPlots/byClass"
# class_pattern <- "tss"
# tssType <- c("Class1","Class2","Class3","Endo")
# dd <-3000
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)
# class_pattern <- "enh"
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)
# class_pattern <- "DHS"
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

# tfType <- "HDAC2"
#  tssType <- c("Class1","Class2","Class3","Endo")
# dd <- c(3000,10000)
# class_pattern <- c("tss","enh","DHS")
# batchOverLapWithOtherFeatures(dd, input.bed.dir, input.bw.path, output.file.dir, class_pattern, tssType, tfType)

# input.bw.path <- "~/Aimin/project/umw_nathan_lawson/Aimin/Alignment"
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/SOX17"
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/SOX17_unsort"

# input.bed.dir <- "/Users/aiminyan/Aimin/project/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/bedFilesforDensityPlots/byClass"
# tfType <- "SOX17"
# tssType <- c("Class1","Class2","Class3","Endo")
# dd <- c(3000,10000)
# class_pattern <- c("tss","enh","DHS")
# batchOverLapWithOtherFeatures(dd, input.bed.dir, input.bw.path, output.file.dir, class_pattern, tssType, tfType)
# tfType <- "GFP"
# batchOverLapWithOtherFeatures(dd, input.bed.dir, input.bw.path, output.file.dir, class_pattern, tssType, tfType)

# tfType <- "SOX17"
# input.bed.dir <- "~/Aimin/project/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/bedFilesforDensityPlots/AllElementsByCellType"
# class_pattern <- c("TSS","Enh")
# tssType <- c("Artery","Common","Vein")
# dd <- c(3000,10000)
# batchOverLapWithOtherFeatures(dd, input.bed.dir, input.bw.path, output.file.dir, class_pattern, tssType, tfType)
# tfType <- "GFP"
# batchOverLapWithOtherFeatures(dd, input.bed.dir, input.bw.path, output.file.dir, class_pattern, tssType, tfType)

batchOverLapWithOtherFeatures <- function(dd, input.bed.dir, input.bw.path, output.file.dir, class_pattern, tssType, tfType) {
  
  null <- lapply(1:length(dd), function(u,input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType){
    
    x <- dd[u]
    null <- lapply(1:length(class_pattern),function(v,input.bed.dir,input.bw.path,output.file.dir,x,class_pattern,tssType,tfTyp){
      
      y <- class_pattern[v]
      
      null <-overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,x,y,tssType,tfType)
      
    },input.bed.dir,input.bw.path,output.file.dir,x,class_pattern,tssType,tfType)
    
  },input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)
  
}

overLapWithOtherFeatures <- function(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType) {
  
  # dd <- 3000
  #  class_pattern <- "tss"
  #  tssType <- "Class1"
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  file.name.4 <- list.files(input.bed.dir,pattern=".bed$",all.files = TRUE,full.names = TRUE,recursive = FALSE,include.dirs = TRUE)
  
  data <- sapply(file.name.4, toGRanges, format="BED")
  
  names(data) <- gsub(".bed", "", basename(names(data)))
  
  #tss <- c("class1tss_2KplMin","class2tss_2KplMin","class3tss_2KplMin","EndoRefTSS_2KplMin")
  
  #tss <- names(data)[grep("tss",names(data),ignore.case = T)]
  
  tss <- names(data)[grep(class_pattern,names(data),ignore.case = T)]
  
  tss.bed <- data[which(names(data) %in% tss)]
  
  feature.center <- lapply(1:length(tss.bed), function(u,tss.bed,dd){
    x <- tss.bed[[u]]
    x
  },tss.bed,dd)
  
  names(feature.center) <- names(tss.bed)
  
  # get the regions of around feature middle +/- dd 
  feature.recentered.l <- lapply(1:length(tss.bed), function(u,tss.bed,dd){
    
    features <- tss.bed[[u]]
    
    wid <- width(features)
    feature.recentered <- feature.center <- features
    
    start(feature.center) <- start(features) + floor(wid/2)
    
    width(feature.center) <- 1
    
    start(feature.recentered) <- start(feature.center) - dd
    end(feature.recentered) <- end(feature.center) + dd
    
    feature.recentered
    
  },tss.bed,dd)
  
  names(feature.recentered.l) <- names(tss.bed)
  
  # feature.test <- tss.bed[[1]]
  # wid.test <- width(feature.test)
  # feature.recentered.test <- feature.center.test <- feature.test
  # start(feature.center.test) <- start(feature.test) + floor(wid.test/2)
  # width(feature.center.test) <- 1
  # start(feature.recentered.test) <- start(feature.center.test) - dd
  # end(feature.recentered.test) <- end(feature.center.test) + dd
  
  cvglists.l <- lapply(1:length(feature.recentered.l), function(u,feature.recentered.l,input.bw.path){
    
    files <- dir(input.bw.path, "bw")
    if(.Platform$OS.type != "windows"){
      cvglists <- sapply(file.path(input.bw.path,files), import, 
                         format="BigWig", 
                         which=feature.recentered.l[[u]], 
                         as="RleList")
    }else{## rtracklayer can not import bigWig files on Windows
      load(file.path(input.bw.path, "cvglist.rds"))
    }
    
  },feature.recentered.l,input.bw.path)
  
  
  grep(tfType,names(unlist(cvglists.l)))
  
  
  #x_name_0 <- tools::file_path_sans_ext(basename(names(cvglists.l[[1]])))  
  #x_name_start <- str_locate_all(x_name_0,"_")[[1]][2,2] + 1
  #x_name_end <- str_locate_all(x_name_0,"_")[[1]][3,2] - 1
  #x_name <- paste0(unique(str_sub(x_name_0,x_name_start,x_name_end)),"_",dd,"_",paste(tss,collapse = "_"))
  
  x_name <- paste0(tfType,"_",dd,"_",paste(tss,collapse = "_"))
  
  saveRDS(cvglists.l,file = file.path(output.file.dir,paste0(x_name,"_cvglist.rds")))
  
  
  
  #lapply(1:length(tss.bed),function(u){
  
  #Class1 <- (cvglists.ul[[1]]+cvglists.ul[[2]])/2
  
  #Class2 <- (cvglists.ul[[3]]+cvglists.ul[[4]])/2
  
  #Class3 <- (cvglists.ul[[5]]+cvglists.ul[[6]])/2
  
  #Endo <- (cvglists.ul[[7]]+cvglists.ul[[8]])/2
  
  #}
  
  #cvglists.ave <- list(Class1=Class1,Class2=Class2,Class3=Class3, Endo=Endo)
  
  cvglists.l <- lapply(1:length(cvglists.l), function(u,cvglists.l){
    
    z <- cvglists.l[[u]][grep(tfType,names(cvglists.l[[u]]))]
    z
  },cvglists.l)
  
  #names(cvglists.l[[4]])
  
  cvglists.ul <- unlist(cvglists.l)
  
  #cvglists.ul <- cvglists.ul[grep(tfType,names(cvglists.ul))]
  
  #grep(tfType,names(unlist(cvglists.l)))
  
  n <- length(cvglists.ul)/length(tss.bed)
  
  if(n==2){
    cvglists.ave <- lapply(1:length(tss.bed),function(u,cvglists.ul){
      
      Class1 <- (cvglists.ul[[2*u-1]]+cvglists.ul[[2*u]])/n
      
      Class1
      
    },cvglists.ul)
  }
  
  if(n==3){
    cvglists.ave <- lapply(1:length(tss.bed),function(u,cvglists.ul){
      
      Class1 <- (cvglists.ul[[3*u-2]]+cvglists.ul[[3*u-1]]+cvglists.ul[[3*u]])/n
      
      Class1
      
    },cvglists.ul)
  }
  
  #cvglists.ave <- list(Class1=Class1,Class2=Class2,Class3=Class3, Endo=Endo)
  
  names(cvglists.ave) <- names(tss.bed)
  
  
  
  
  # sig1 <- featureAlignedSignal(cvglists.ave, feature.center[[1]], upstream=dd, downstream=dd) 
  # sigs.log2 <- lapply(sig1, function(.ele) log2(.ele+1))   
  # out.1 <- featureAlignedDistribution(sigs.log2, 
  #                                   feature.center[[1]],upstream=dd, downstream=dd,
  #                                   zeroAt=.5,type="l", 
  #                                   ylab="Averaged coverage")
  # 
  # 
  # sig1 <- featureAlignedSignal(cvglists.ave, feature.center[[2]], upstream=dd, downstream=dd) 
  # sigs.log2 <- lapply(sig1, function(.ele) log2(.ele+1))
  # out.2 <- featureAlignedDistribution(sigs.log2, 
  #                                   feature.center[[2]],upstream=dd, downstream=dd,
  #                                   zeroAt=.5,type="l", 
  #                                   ylab="Averaged coverage",add=TRUE)
  # 
  # sig1 <- featureAlignedSignal(cvglists.ave, feature.center[[3]], upstream=dd, downstream=dd) 
  # sigs.log2 <- lapply(sig1, function(.ele) log2(.ele+1))
  # out.3 <- featureAlignedDistribution(sigs.log2, 
  #                                   feature.center[[3]],upstream=dd, downstream=dd,
  #                                   zeroAt=.5,type="l", 
  #                                   ylab="Averaged coverage",add=TRUE)
  # 
  # 
  # sig1 <- featureAlignedSignal(cvglists.ave, feature.center[[4]], upstream=dd, downstream=dd) 
  # sigs.log2 <- lapply(sig1, function(.ele) log2(.ele+1))
  # out.4 <- featureAlignedDistribution(sigs.log2, 
  #                                   feature.center[[4]],upstream=dd, downstream=dd,
  #                                   zeroAt=.5,type="l", 
  #                                   ylab="Averaged coverage",add=TRUE)
  # 
  # out <- cbind(out.1$density[,1],out.2$density[,2],out.3$density[,3],out.4$density[,4])
  
  XX <- lapply(1:length(feature.center),function(u,feature.center,cvglists.ave,dd){
    
    # u <- 1
    sig1 <- featureAlignedSignal(cvglists.ave, feature.center[[u]], upstream=dd, downstream=dd) 
    sigs.log2 <- lapply(sig1, function(.ele) log2(.ele+1))
    out.4 <- featureAlignedDistribution(sigs.log2, 
                                        feature.center[[u]],upstream=dd, downstream=dd,
                                        zeroAt=.5,type="l", 
                                        ylab="Averaged coverage",add=TRUE)
    
    x <- list(density=out.4$density[,u],grWidAt=out.4$grWidAt,grWidLab=out.4$grWidLab)
    x 
  },feature.center,cvglists.ave,dd)
  
  den.l <-lapply(1:length(XX),function(u,XX){
    Z <- XX[[u]]$density
    Z
  },XX)
  out <- do.call(cbind,den.l)
  
  grWidAt.l <-lapply(1:length(XX),function(u,XX){
    Z <- XX[[u]]$grWidAt
    Z
  },XX)
  grWidAt <- do.call(cbind,grWidAt.l)
  
  grWidLab.l <-lapply(1:length(XX),function(u,XX){
    Z <- XX[[u]]$grWidLab
    Z
  },XX)
  grWidLab <- do.call(cbind,grWidLab.l)
  
  
  head(out)
  colnames(out) <- tssType
  
  drowDensity <- function(density,grWidAt,grWidLab,...) {
    dots <- list(...)
    matplot(density, ..., xaxt="n")
    axis(1, at = grWidAt, labels = grWidLab)
    lty <- if(!is.null(dots$lty)) dots$lty else 1:5
    lwd <- if(!is.null(dots$lwd)) dots$lwd else 1
    col <- if(!is.null(dots$col)) dots$col else 1:6
    legend("topright", legend=colnames(density), col=col,
           lty=lty, lwd=lwd)
  }
  
  jpeg(file.path(output.file.dir,paste0(x_name,"_density.jpeg")))
  drowDensity(out,grWidAt[,1],grWidLab[,1],type="l", lty=1,
              ylab="log2(averaged coverage)",xlab=paste0("Distance to ",class_pattern,"(bp)"),main=tfType)
  dev.off()
  
  x.pos <- seq(-dd,dd,length.out = dim(out)[1])
  out.m <- data.frame(x.pos=x.pos,out)
  
  outLong <- melt(data          = out.m,
                  id.vars       = c("x.pos"),
                  measure.vars  = tssType,
                  variable.name = "TSSByClass",
                  value.name    = "value")
  
  ggplot(outLong, aes(x.pos, value,colour=TSSByClass)) + 
    geom_line() + labs(x="Distance to TSS(bp)",y="log2(averaged coverage)",title=tfType) +
    theme(plot.title=element_text(hjust=0.5)) + scale_x_continuous(name=paste0("Distance to ",class_pattern,"(bp)"), limits=c(-dd, dd))
  
  write.table(out.m,file = file.path(output.file.dir,paste0(x_name,"_Density_matrix_short",".txt")),
              append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T)
  
  write.table(outLong,file = file.path(output.file.dir,paste0(x_name,"_Density_matrix_long",".txt")),
              append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T)
  
  null <- lapply(1:length(cvglists.l), function(u,cvglists.l,feature.center,dd,n,tfType){
    #u <- 1
    sig1.4.heatmap <- featureAlignedSignal(cvglists.l[[u]], feature.center[[u]], upstream=dd, downstream=dd) 
    sig1.4.heatmap.log2 <- lapply(sig1.4.heatmap, function(.ele) log2(.ele+1))
    names(sig1.4.heatmap.log2) <- paste0(tfType,"_",seq(1,n))
    
    #u <- 1
    #x_name_0 <- tools::file_path_sans_ext(basename(names(cvglists.l[[u]])))  
    #x_name_start <- str_locate_all(x_name_0,"_")[[1]][2,2] + 1
    #x_name_end <- str_locate_all(x_name_0,"_")[[1]][3,2] - 1
    
    x_name <- paste0(tfType,"_",names(feature.center)[u])
    
    jpeg(file.path(output.file.dir,paste0(x_name,"_",dd,"_heatmap.jpeg")))
    featureAlignedHeatmap(sig1.4.heatmap.log2, feature.center[[u]],sortBy=NULL,upstream=dd, downstream=dd,zeroAt=.5,res=300)
    dev.off()
    
  },cvglists.l,feature.center,dd,n,tfType)
  
}

# path <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/H3K27me3_bigwig"
# path <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/NR2F2kd_bigwig"

getCvg <- function(path,output.file.dir) {
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  # K27me3.bw.dir <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/H3K27me3_bigwig"
  # NR2F2kd.bw.dir <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/NR2F2kd_bigwig"
  
  files <- dir(path, "bw")
  
  if(.Platform$OS.type != "windows"){
    cvglists <- sapply(file.path(path, files), import, 
                       format="BigWig", 
                       as="RleList")
  }else{## rtracklayer can not import bigWig files on Windows
    load(file.path(path, "cvglist.rds"))
  }
  
  saveRDS(cvglists,file = file.path(output.file.dir,"cvglist.rds"))
  
}

# nathan_lawson new data

# input.bed.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/newBedFiles"
# input.bw.path <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/NR2F2kd_bigwig"
# dd <- 3000
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/NR2F2_new"
# null <- overLapWithOtherFeatures1(input.bed.dir,input.bw.path,output.file.dir,dd)

overLapWithOtherFeatures1 <- function(input.bed.dir,input.bw.path,output.file.dir,dd) {
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  file.name.4 <- list.files(input.bed.dir,pattern=".bed$",all.files = TRUE,full.names = TRUE,recursive = FALSE,include.dirs = TRUE)
  data <- sapply(file.name.4, toGRanges, format="BED")
  
  names(data) <- gsub(".bed", "", basename(names(data)))
  
  generateMatrix4heatmap <- function(class_pattern,data,dd,input.bw.path,sortBy) {
    
    tss <- names(data)[grep(class_pattern,names(data),ignore.case = T)]
    print(tss)
    
    tss.bed <- data[which(names(data) %in% tss)]
    
    feature.center <- lapply(1:length(tss.bed), function(u,tss.bed,dd){
      x <- tss.bed[[u]]
      x
    },tss.bed,dd)
    
    names(feature.center) <- names(tss.bed)
    
    # get the regions of around feature middle +/- dd 
    feature.recentered.l <- lapply(1:length(tss.bed), function(u,tss.bed,dd){
      
      features <- tss.bed[[u]]
      
      wid <- width(features)
      feature.recentered <- feature.center <- features
      
      start(feature.center) <- start(features) + floor(wid/2)
      
      width(feature.center) <- 1
      
      start(feature.recentered) <- start(feature.center) - dd
      end(feature.recentered) <- end(feature.center) + dd
      
      feature.recentered
      
    },tss.bed,dd)
    names(feature.recentered.l) <- names(tss.bed)
    
    cvglists.l <- lapply(1:length(feature.recentered.l), function(u,feature.recentered.l,input.bw.path){
      
      files <- dir(input.bw.path, "bw")
      if(.Platform$OS.type != "windows"){
        cvglists <- sapply(file.path(input.bw.path,files), import, 
                           format="BigWig", 
                           which=feature.recentered.l[[u]], 
                           as="RleList")
      }else{## rtracklayer can not import bigWig files on Windows
        load(file.path(input.bw.path, "cvglist.rds"))
      }
      
    },feature.recentered.l,input.bw.path)
    
    tfType <- "TRC-NS"
    cvglists.ck <- lapply(1:length(cvglists.l), function(u,cvglists.l){
      z <- cvglists.l[[u]][grep(tfType,names(cvglists.l[[u]]))]
      z
    },cvglists.l)
    cvglists.ul.ck <- unlist(cvglists.ck)
    
    
    print(names(cvglists.ul.ck))
    
    n <- length(cvglists.ul.ck)
    
    sig1.4.heatmap.ck <- featureAlignedSignal(cvglists.ck[[1]], feature.center[[1]], upstream=dd, downstream=dd)
    sig1.4.heatmap.log2.ck <- lapply(sig1.4.heatmap.ck, function(.ele) log2(.ele+1))
    names(sig1.4.heatmap.log2.ck) <- paste0(tfType,"_",seq(1,n))
    print(names(sig1.4.heatmap.log2.ck))
    
    tfType <- "shNR2F2"
    cvglists.tr <- lapply(1:length(cvglists.l), function(u,cvglists.l){
      z <- cvglists.l[[u]][grep(tfType,names(cvglists.l[[u]]))]
      z
    },cvglists.l)
    cvglists.ul.tr <- unlist(cvglists.tr)
    print(names(cvglists.ul.tr))
    n <- length(cvglists.ul.tr)
    
    sig1.4.heatmap.tr <- featureAlignedSignal(cvglists.tr[[1]], feature.center[[1]], upstream=dd, downstream=dd)
    sig1.4.heatmap.log2.tr <- lapply(sig1.4.heatmap.tr, function(.ele) log2(.ele+1))
    names(sig1.4.heatmap.log2.tr) <- paste0(tfType,"_",seq(1,n))
    print(names(sig1.4.heatmap.log2.tr))
    
    dim(sig1.4.heatmap.log2.tr[[1]])
    
    #sortBy = "Trt"
    sig1.4.heatmap.log2.ck.tr.l <- lapply(1:length(sig1.4.heatmap.log2.tr),function(u,sig1.4.heatmap.log2.ck,sig1.4.heatmap.log2.tr,sortBy){
      if(sortBy=="Trt"){
        density.sum <- apply(sig1.4.heatmap.log2.tr[[u]],1,sum)
        sig1.4.heatmap.log2.tr.sorted<- sig1.4.heatmap.log2.tr[[u]][order(density.sum,decreasing = T),]        
        sig1.4.heatmap.log2.ck.sorted <- sig1.4.heatmap.log2.ck[[u]][match(rownames(sig1.4.heatmap.log2.tr.sorted),rownames(sig1.4.heatmap.log2.ck[[u]])),]
        sig1.4.heatmap.log2.ck.tr <- cbind(sig1.4.heatmap.log2.ck.sorted,sig1.4.heatmap.log2.tr.sorted)
      }else{
        density.sum <- apply(sig1.4.heatmap.log2.ck[[u]],1,sum)
        sig1.4.heatmap.log2.ck.sorted<- sig1.4.heatmap.log2.ck[[u]][order(density.sum,decreasing = T),]        
        sig1.4.heatmap.log2.tr.sorted <- sig1.4.heatmap.log2.tr[[u]][match(rownames(sig1.4.heatmap.log2.ck.sorted),rownames(sig1.4.heatmap.log2.tr[[u]])),]
        sig1.4.heatmap.log2.ck.tr <- cbind(sig1.4.heatmap.log2.ck.sorted,sig1.4.heatmap.log2.tr.sorted)
      }
      sig1.4.heatmap.log2.ck.tr
    },sig1.4.heatmap.log2.ck,sig1.4.heatmap.log2.tr,sortBy)
    
    sig1.4.heatmap.log2.ck.tr.l
    
  }
  
  heatmap4UpNoChangeDown <- function(class_patterns,file_name,data, dd, input.bw.path, output.file.dir) {
    class_pattern <- class_patterns[1]
    sig1.4.heatmap.log2.ck.tr.no.change <- generateMatrix4heatmap(class_pattern,data,dd,input.bw.path,sortBy="Trt") 
    
    class_pattern <- class_patterns[2]
    sig1.4.heatmap.log2.ck.tr.Up <- generateMatrix4heatmap(class_pattern,data,dd,input.bw.path,sortBy="Trt")
    
    class_pattern <- class_patterns[3]
    sig1.4.heatmap.log2.ck.tr.Down <- generateMatrix4heatmap(class_pattern,data,dd,input.bw.path,sortBy="CK")
    
    ht_list <- lapply(1:length(sig1.4.heatmap.log2.ck.tr.no.change), function(u,sig1.4.heatmap.log2.ck.tr.no.change,sig1.4.heatmap.log2.ck.tr.Up,sig1.4.heatmap.log2.ck.tr.Down){
      
      sig1.4.heatmap.log2.ck.tr.1 <- rbind(sig1.4.heatmap.log2.ck.tr.no.change[[u]],sig1.4.heatmap.log2.ck.tr.Up[[u]],sig1.4.heatmap.log2.ck.tr.Down[[u]])
      
      z.mat.0 <- t(scale(t(sig1.4.heatmap.log2.ck.tr.1), center=TRUE, scale=TRUE))
      
      n.no.change <- dim(sig1.4.heatmap.log2.ck.tr.no.change[[u]])[1]
      n.Up <- dim(sig1.4.heatmap.log2.ck.tr.Up[[u]])[1]
      n.Down <- dim(sig1.4.heatmap.log2.ck.tr.Down[[u]])[1]
      
      ht <- Heatmap(z.mat.0, name = paste0("z-score","-rep",u),
                    col = rev(redgreen(30))[-seq(35, 35)],            
                    show_row_name = FALSE,
                    cluster_columns = F,
                    cluster_rows = F,column_names_gp = gpar(fontsize = 6),column_split = rep(c("Control", "NR2F2KD"),c(100,100)),row_split = factor(rep(c("NoChange","Up","Down"),c(n.no.change,n.Up,n.Down)),levels = c("Up","NoChange","Down")))
      
      ht
      
    },sig1.4.heatmap.log2.ck.tr.no.change,sig1.4.heatmap.log2.ck.tr.Up,sig1.4.heatmap.log2.ck.tr.Down)
    
    htList = ht_list[[1]] + ht_list[[2]] + ht_list[[3]]
    
    gb = grid.grabExpr(draw(htList))
    g <- plot_grid(gb)
    save_plot(file.path(output.file.dir,paste0(file_name,"_heatmap.png")),g,base_width = 10,base_height = 40)
  }
  
  class_patterns <- c("Com_NoChange","Com_K27acUp","Com_K27acDwn")
  file_name <- "Com_3_reps"
  heatmap4UpNoChangeDown(class_patterns,file_name,data, dd, input.bw.path, output.file.dir)
  
  class_patterns <- c("Art_NoChange","Art_K27acUp","Art_K27acDown")
  file_name <- "Art_3_reps"
  heatmap4UpNoChangeDown(class_patterns,file_name,data, dd, input.bw.path, output.file.dir)
  
  class_patterns <- c("Ve_NoChange","Ve_K27acUp","Ve_K27acDwn")
  file_name <- "Ve_3_reps"
  heatmap4UpNoChangeDown(class_patterns,file_name,data, dd, input.bw.path, output.file.dir)
  
}

featureAlignedHeatmap <- function(cvglists, feature.gr, upstream, downstream, 
                                  zeroAt, n.tile=100,
                                  annoMcols=c(), sortBy=names(cvglists)[1],
                                  color=colorRampPalette(c("yellow", "red"))(50),
                                  lower.extreme, upper.extreme,
                                  margin=c(0.1, 0.01, 0.15, 0.1), gap=0.01, 
                                  newpage=TRUE, gp=gpar(fontsize=10),
                                  ...){
  
  
  cat("local one \n")
  
  if(!missing(lower.extreme)){
    stopifnot(is.numeric(lower.extreme))
  }
  if(!missing(upper.extreme)){
    stopifnot(is.numeric(upper.extreme))
  }
  stopifnot(is(gp, "gpar"))
  stopifnot(is.logical(newpage))
  stopifnot(is(feature.gr, "GRanges"))
  
  #grWidr <- unique(width(feature.gr))
  grWidr <- width(feature.gr)
  if(missing(upstream) || missing(downstream)){
    if(length(grWidr)!=1){
      stop("The width of feature.gr is not identical.")
    }
    if(missing(zeroAt)) {
      zeroAt <- 0.5
      message("zero is set as the center of the feature.gr")
    }
  }else{
    if(!is.numeric(upstream) || !is.numeric(downstream)){
      stop("upstream and downstream must be integers")
    }
    if(upstream<0 || downstream<0){
      stop("upstream and downstream must be not less than 0")
    }
    upstream <- as.integer(upstream)
    downstream <- as.integer(downstream)
    if(length(grWidr)!=1){
      start(feature.gr) <- start(feature.gr)+floor(grWidr/2)
      width(feature.gr) <- 1
      warning("feature.gr is set to the center of feature.gr")
      # print(feature.gr)
    }
    zeroAt <- upstream/(upstream + downstream)
    feature.gr <- promoters(feature.gr, upstream = upstream,
                            downstream = downstream + 1)
    #print(feature.gr)
    grWidr <- unique(width(feature.gr))
    #grWidr <- width(feature.gr)
  }
  stopifnot(is.numeric(zeroAt))
  stopifnot(zeroAt>=0)
  if(zeroAt<=1){
    
    cat(grWidr," ", zeroAt,"\n")
    
    zero <- round(grWidr*zeroAt)
    print(zero)
  }else{
    zero <- round(zeroAt)
  }
  
  grWid <- c(0, grWidr) - zero
  print(grWid)
  grWidLab <- grid.pretty(grWid)
  print(grWidLab)
  grWidAt <- (grWidLab+zero)/grWidr
  print(grWidAt)
  if(inherits(cvglists, c("SimpleRleList", "RleList", "CompressedRleList"))){
    
    cat("Test this \n")
    cvglistsName <- substitute(deparse(cvglists))
    cvglists <- list(cvglists)
    names(cvglists) <- cvglistsName
  }
  if(!is(cvglists, "list")){
    stop("cvglists must be output of featureAlignedSignal or 
         a list of SimpleRleList or RleList")
  }
  cls <- sapply(cvglists, class)
  if(all(cls=="matrix")){
    cov <- cvglists
    print(cov)
    if(ncol(cov[[1]])!=n.tile){
      stop("n.tile must keep same as featureAligendSignal.")
    }
  }else{
    if(any(!cls %in% c("SimpleRleList", "RleList", "CompressedRleList")))
      stop("cvglists must be a list of SimpleRleList or RleList")
    cov <- featureAlignedSignal(cvglists, feature.gr, n.tile=n.tile)
  }
  
  print(feature.gr)
  
  annoMcols <- annoMcols[annoMcols %in% colnames(mcols(feature.gr))]
  if(length(sortBy)>0){
    covTotalCoverage <- 
      do.call(cbind, lapply(cov, rowSums, na.rm=TRUE))
    
    print(covTotalCoverage)
    
    colnames(covTotalCoverage) <- names(cvglists)
    if(length(annoMcols)>0){
      sortMatrix <- mcols(feature.gr)[, annoMcols, drop=FALSE]
      for(i in seq.int(ncol(sortMatrix))){
        if(!is.numeric(sortMatrix[,i])){
          sortMatrix[,i] <- 
            factor(as.character(sortMatrix[,i]), 
                   levels=rev(unique(as.character(sortMatrix[, i]))))
        }
      }
      if(all(sortBy %in% colnames(mcols(feature.gr)))){
        sortMatrix <- 
          cbind(sortMatrix, 
                as.data.frame(mcols(feature.gr)[, sortBy, drop=FALSE]))
      }else{
        if(all(sortBy %in% colnames(covTotalCoverage))){
          sortMatrix <- 
            cbind(sortMatrix,
                  data.frame(covTotalCoverage[, sortBy, drop=FALSE]))
        }else{
          stop("sortBy is incorrect.")
        }
      }
    }else{
      if(all(sortBy %in% colnames(mcols(feature.gr)))){
        sortMatrix <- 
          as.data.frame(mcols(feature.gr)[, sortBy, drop=FALSE])
      }else{
        if(all(sortBy %in% colnames(covTotalCoverage))){
          
          sortMatrix <- data.frame(covTotalCoverage[, sortBy, drop=FALSE])
          
          print(sortMatrix)
          
        }else{
          stop("sortBy is incorrect.")
        }
      }
    }
    names(sortMatrix) <- paste("V", 1:ncol(sortMatrix))
    soid <- do.call(order, c(as.list(sortMatrix), decreasing = TRUE))
    cov <- lapply(cov, function(.ele) .ele[soid, ])
    feature.gr <- feature.gr[soid]
  }
  
  #covert to color
  if(missing(lower.extreme)){
    lower.extreme <- rep(0, length(cvglists))
  }
  if(missing(upper.extreme)){
    lim=list()
  }else{
    upper.extreme <- 
      rep(upper.extreme, length(cvglists))[1:length(cvglists)]
    lim=as.list(data.frame(rbind(lower.extreme, upper.extreme)))
    names(lim) <- NULL
  }
  
  if(length(lim)==0){
    lim <- sapply(cov, range, na.rm=TRUE, simplify=FALSE)
  }
  if(length(lim)==1){
    lim <- rep(lim, length(cov))
    legend <- 1
  }else{
    legend <- length(cov)
  }
  if(any(sapply(lim, length)<2)){
    stop("each limit should be a numeric vector with length 2")
  }
  cov <- mapply(function(.ele, .lim){
    if(min(.lim)==max(.lim)){
      .lim <- seq(min(.lim), max(.lim)+1, length.out=length(color)+1)
    }else{
      .lim <- seq(min(.lim), max(.lim), length.out=length(color)+1)
    }
    .lim[1] <- min(.lim[1], 0)
    .lim[length(.lim)] <- max(max(.lim), .Machine$integer.max)
    t(apply(.ele, 1, 
            cut, breaks=.lim, labels=color, include.lowest=TRUE))
  }, cov, lim, SIMPLIFY=FALSE)
  # library(grid)
  # grid.raster
  if(newpage) grid.newpage()
  allGrob <- NULL
  height=1-margin[3]-margin[1]
  width=1-margin[4]-margin[2]
  ## draw y labels
  isFloat <- function(x){
    if(is.numeric(x)){
      return(x!=floor(x))
    }
    return(FALSE)
  }
  areColors <- function(x) {
    sapply(x, function(X) {
      if(isFloat(X)) return(FALSE)
      if(is.numeric(X)){
        if(X<9 && X>0){
          return(TRUE)
        }else{
          return(FALSE)
        }
      }
      tryCatch(is.matrix(col2rgb(X)), 
               error = function(e) FALSE)
    })
  }
  colorGroup <- 
    list(c("#295AA5", "#52BDBD", "#FFAD6B", "#F784A5", "#CEB55A", 
           "#31B5D6", "#BD2119", "#F79463", "#FFE608", "#D69419", 
           "#29B5CE", "#737BB5"), 
         c("#4A4221", "#A5BD42", "#8CC584", "#CEE6D6", "#C5D64A", 
           "#FFF79C", "#84AD31", "#DBDECE", "#D6AD4A", "#4A4229", 
           "#637B7B", "#525229", "#F7F7AD", "#EFAD63"), 
         c("#6B5252", "#DE6B73", "#000000", "#4A4221", "#3AA55A", 
           "#DEE6C5", "#102119", "#DEE6CE", "#738CC5", "#31317B", 
           "#F7F7AD", "#8C9CA5", "#732163", "#001000"), 
         c("#192919", "#8C946B", "#6B1921", "#F7F7AD", "#4A4229", 
           "#CEA54A", "#000000", "#6B736B", "#524A29", "#423119", 
           "#213129", "#7B4229", "#7B2121", "#C57B4A"), 
         c("#9D2932", "#EFDEB0", "#789262", "#494166", "#A88462", 
           "#D9B612", "#177CB0", "#F3F8F1", "#424B50", "#B36D61", 
           "#C3272B", "#A98175", "#D4F2E8", "#FF3300", "#76664D", 
           "#ED5736", "#A3E2C5", "#415065", "#D7ECF1"), 
         c("#009C8C", "#0094AD", "#BDDEC5", "#00ADCE", "#DEFF19", 
           "#D6EFB5", "#EFF7E6", "#42CE73", "#42A5CE", "#BDFFF7", 
           "#DEFFDE", "#4A6BCE"), 
         c("#D6B5BD", "#EFD6BD", "#FFE6D6", "#DEE6CE", "#F7F7C5", 
           "#D6E6B5", "#73BD8C", "#BDDEDE", "#DEBD8C", "#E6A5BD", 
           "#F7ADAD", "#FFE6D6", "#F7F7D6", "#E6BDBD"), 
         c("#E6E6BD", "#637342", "#425A52", "#298473", "#BDDEC5", 
           "#527384", "#73847B", "#CED6AD", "#94D6CE", "#527384", 
           "#4A6B5A", "#318442", "#C5DEDE"), 
         c("#F7846B", "#F79C9C", "#CEE6E6", "#FFEFBD", "#EFEFBD", 
           "#FFE6DE", "#FFBDA5", "#FFEF00", "#63BD84", "#D6E6AD", 
           "#FFF7D6", "#FFFFEF", "#FFE6DE", "#FFBD9C", "#EF9C94"), 
         c("#4A4229", "#212108", "#6B8484", "#EFE6BD", "#E68C3A", 
           "#A56352", "#6B2921", "#525229", "#844229", "#8C524A", 
           "#E68C3A", "#D6AD4A", "#212108", "#524A3A"), 
         c("#312184", "#AD2919", "#EF3129", "#D69C4A", "#FF6B00", 
           "#6B4A31", "#CE3121", "#AD2919", "#5A3A5A", "#292931", 
           "#3A3184", "#9C5229", "#EFCE19"), 
         c("#6373B5", "#FF5200", "#4AF54A", "#FFCE10", "#EF3A21", 
           "#EF4A94", "#5A63AD", "#FFEDA5", "#EF3121", "#FFC510", 
           "#FF5200"), 
         c("#FFB6C1", "#DC143C", "#DB7093", "#FF69B4", "#FF1493", 
           "#C71585", "#DA70D6", "#FF00FF", "#8B008B", "#9400D3", 
           "#4B0082", "#FFA07A", "#FF7F50", "#FF4500", "#E9967A", 
           "#FF6347", "#F08080", "#CD5C5C", "#FF0000", "#A52A2A", 
           "#B22222", "#8B0000"), 
         c("#7B68EE", "#6A5ACD", "#483D8B", "#E6E6FA", "#0000FF", 
           "#191970", "#4169E1", "#6495ED", "#B0C4DE", "#1E90FF", 
           "#87CEFA", "#00BFFF", "#00FFFF"), 
         c("#008B8B", "#48D1CC", "#7FFFD4", "#66CDAA", "#00FA9A", 
           "#00FF7F", "#3CB371", "#2E8B57", "#90EE90", "#32CD32", 
           "#00FF00", "#008000", "#006400", "#ADFF2F", "#556B2F", 
           "#9ACD32", "#6B8E23"), 
         c("#FAFAD2", "#FFFF00", "#808000", "#BDB76B", "#FFFACD", 
           "#EEE8AA", "#F0E68C", "#FFD700", "#DAA520", "#B8860B", 
           "#F5DEB3", "#FFA500", "#FFEBCD", "#FFDEAD", "#D2B48C", 
           "#FF8C00", "#CD853F", "#8B4513"),
         c("#F5F5F5", "#DCDCDC", "#D3D3D3", "#C0C0C0", "#A9A9A9", 
           "#808080", "#696969", "#5C5C5C", "#4D4D4D", "#333333", 
           "#1A1A1A", "#000000"))
  if(length(annoMcols)>0){
    ci <- 1
    ytop <- 0.75-margin[3]
    mcols <- mcols(feature.gr)
    vp.annoMcols <- viewport(x=.5, y=height/2+margin[1], 
                             width=1, height=height,
                             name="vp.annoMcols")
    for(i in annoMcols){
      mc <- as.character(mcols[, i])
      mc.name <- ifelse(is.numeric(i), colnames(mcols)[i], i)
      if(!all(areColors(mcols[, i]))){
        if(is.numeric(mcols[, i])){
          mc.range <- range(mcols[, i])
          mc.color <- range(colorGroup[[ci]])
          mc.color <- colorRampPalette(mc.color)(100)
          mc <- cut(mcols[, i], breaks = 100, labels = mc.color)
          mc <- as.character(mc)
          mc.label <- sort(grid.pretty(mc.range), decreasing = TRUE)
          mc.label.color <-
            mc.color[round(100*(mc.label-min(mc.range))/
                             diff(mc.range))]
          mc.label <- mc.label[!is.na(mc.label.color)]
          mc.label.color <- mc.label.color[!is.na(mc.label.color)]
        }else{
          mc.label <- unique(mc)
          mc <- factor(mc, levels=mc.label)
          mc.label.color <- rep(colorGroup[[ci]], 
                                length(mc.label))[1:length(mc.label)]
          levels(mc) <- mc.label.color
          mc <- as.character(mc)
        }
        
        ci <- ci+1
      }else{
        mc.label <- unique(mcols[, i])
        if(is.numeric(mcols[, i])){
          mc <- col2rgb(mc);
          mc <- rgb(red=mc[1, ],
                    green=mc[2, ],
                    blue=mc[3, ], maxColorValue = 255)
        }
        mc.label.color <- unique(mc)
      }
      ## mc is color
      ## mc.label is legends
      ## plot legned
      ht.annoMcols.legend <- length(mc.label) * 0.03
      wd.annoMcols.legend <- margin[4]
      vp.annoMcols.legend <- 
        viewport(x=1-wd.annoMcols.legend+gap, 
                 y=ytop-ht.annoMcols.legend/2, 
                 width=wd.annoMcols.legend, 
                 height=ht.annoMcols.legend,
                 just=0,
                 name="vp.annoMcols.legend")
      annoMcols.legend <- legendGrob(labels=mc.label, 
                                     ncol=1,
                                     pch=15,
                                     gp=gpar(col=mc.label.color))
      allGrob <- gList(allGrob, 
                       gTree(children=gList(annoMcols.legend),
                             vp=vp.annoMcols.legend,
                             name=paste("gTree.annoMcols.legend",
                                        mc.name,
                                        sep=".")))
      ytop <- ytop-ht.annoMcols.legend-gap
      
      ## add annoMcols
      vp.annoMcols.sub <- viewport(x=width - .0075,
                                   y=.5, 
                                   width=.015,
                                   height=1,
                                   name=paste("vp.annoMcols",
                                              mc.name,
                                              sep="."))
      raster.annoMcols.sub <- 
        rasterGrob(mc,
                   width=1, height=1,
                   interpolate=FALSE,
                   name=paste("raster.annoMcols",
                              mc.name,
                              sep="."),
                   vp=vp.annoMcols.sub)
      xaxis.annoMcols.sub <- 
        xaxisGrob(at=.5,
                  label=mc.name,
                  gp=gp,
                  edits=gEdit("labels", rot=90),
                  name=paste("xaxis.annoMcols",
                             mc.name,
                             sep="."),
                  vp=vp.annoMcols.sub)
      allGrob <- gList(allGrob, 
                       gTree(children=gList(raster.annoMcols.sub,
                                            xaxis.annoMcols.sub),
                             vp=vp.annoMcols,
                             name=paste("gTree.annoMcols",
                                        mc.name,
                                        sep=".")))
      width <- width - .015
    }
    width <- width - gap
  }
  
  vp.heatmap <- viewport(x=width/2+margin[2], y=height/2+margin[1], 
                         width=width, height=height,
                         name="vp.heatmap")
  l <- length(cov)
  wid <- 1/l
  x <- -wid/2
  y <- 1
  for(i in 1:l){
    x <- x+wid
    x.name <- names(cvglists)[i]
    vp.heatmap.sub <- viewport(x=x, width=wid-gap, height=y,
                               name=paste("vp.heatmap",
                                          x.name,
                                          sep="."))
    raster.heatmap.sub <- rasterGrob(cov[[i]], 
                                     width=1, height=1, 
                                     interpolate=FALSE,
                                     name=paste("raster.heatmap",
                                                x.name,
                                                sep="."),
                                     vp=vp.heatmap.sub)
    
    grWidLab.this <- grWidLab
    if(i==1){
      if(i!=l){
        if(length(grWidAt)>1){
          grWidLab.this[length(grWidLab)] <- " "
        }
      }
    }else{
      if(i==l){##i!=1
        if(length(grWidAt)>1){
          grWidLab.this[1] <- " "
        }
      }else{## middle
        if(length(grWidAt)>2){
          grWidLab.this[1] <- " "
          grWidLab.this[length(grWidLab)] <- " "
        }
      }
    }
    xaxis.heatmap.sub <- xaxisGrob(at=grWidAt,
                                   label=grWidLab.this,
                                   gp=gp,
                                   edits=gEdit("labels", rot=90),
                                   name=paste("xaxis.heatmap",
                                              x.name,
                                              sep="."),
                                   vp=vp.heatmap.sub)
    allGrob <- gList(allGrob, 
                     gTree(children=gList(raster.heatmap.sub,
                                          xaxis.heatmap.sub),
                           vp=vp.heatmap,
                           name=paste("gTree.heatmap",
                                      x.name,
                                      sep=".")))
  }
  
  ## add heatmap legend
  if(legend==1){
    vp.legend <- viewport(x=1.01-margin[4]+gap, y=.875-margin[3], 
                          width=.02, height=.25,
                          name="vp.legend")
    raster.legend <- rasterGrob(matrix(rev(color), ncol=1), 
                                width=1, height=1,
                                name="raster.legend")
    .lim <- lim[[1]]
    label <- grid.pretty(.lim)
    ll <- length(label)
    if(ll>=2){
      if(ll%%2==0){
        label <- label[c(1, ll)]
      }else{
        label <- label[c(1, ceiling(ll/2), ll)]
      }
      at <- (label-min(.lim))/abs(diff(.lim))
      yaxis.legend <- yaxisGrob(at=at,
                                label=label,
                                main=FALSE,
                                gp=gp,
                                name="yaxis.legend")
    }
    allGrob <- gList(allGrob, 
                     gTree(children=gList(raster.legend,
                                          yaxis.legend),
                           vp=vp.legend,
                           name="gTree.legend"))
    
    ## draw x labels
    vp.xlabels <- viewport(x=width/2+margin[2], 
                           y=1-margin[3]/2+gap, 
                           width=width, height=margin[3]-gap,
                           name="vp.xlabels")
    x <- -wid/2
    h <- .5
    for(i in 1:l){
      x <- x+wid
      x.name <- names(cvglists)[i]
      text.xlabels.sub <- 
        textGrob(label = x.name,
                 gp=gp,
                 just=1, hjust=0,
                 vp=viewport(x=x, y=.5, width=wid, height=1, 
                             name=paste("vp.xlabels", x.name, sep=".")),
                 name=paste("text.xlabels", x.name, sep="."))
      allGrob <- gList(allGrob, 
                       gTree(children=gList(text.xlabels.sub),
                             vp=vp.xlabels,
                             name=paste("gTree.xlabels",
                                        x.name,
                                        sep=".")))
    }
  }else{
    vp.legend <- viewport(x=width/2+margin[2], 
                          y=1-margin[3]/3+gap, 
                          width=width, height=margin[3]*2/3-gap,
                          name="vp.legend")
    
    vp.xlabels <- viewport(x=width/2+margin[2], 
                           y=1-margin[3]*5/6+gap, 
                           width=width, height=margin[3]/3-gap,
                           name="vp.xlabels")
    
    x <- -wid/2
    h <- .5
    for(i in 1:l){
      x <- x+wid
      x.name=names(cvglists)[i]
      vp.legend.sub <- 
        viewport(x=x, y=0.25, width=wid/2, height=0.2,
                 name=paste("vp.legend",
                            x.name,
                            sep="."))
      raster.legend.sub <- 
        rasterGrob(matrix(color, nrow=1), 
                   width=1, height=1,
                   name=paste("raster.legend", x.name, sep="."),
                   vp=vp.legend.sub)
      .lim <- lim[[i]]
      label <- grid.pretty(.lim)
      ll <- length(label)
      xaxis.legned.sub <- NULL
      if(ll>=2){
        if(ll%%2==0){
          label <- label[c(1, ll)]
        }else{
          label <- label[c(1, ceiling(ll/2), ll)]
        }
        at <- (label-min(.lim))/abs(diff(.lim))
        xaxis.legned.sub <- xaxisGrob(at=at,
                                      label=label,
                                      main=FALSE,
                                      gp=gp,
                                      name=paste("xaxis.legend",
                                                 x.name,
                                                 sep="."),
                                      vp=vp.legend.sub)
      }
      allGrob <- gList(allGrob, 
                       gTree(children=gList(raster.legend.sub,
                                            xaxis.legned.sub),
                             vp=vp.legend,
                             name=paste("gTree.legend",
                                        x.name,
                                        sep=".")))
      
      ## draw x labels
      text.xlabels.sub <- 
        textGrob(label = x.name,
                 gp=gp,
                 vp=viewport(x=x, y=.5, width=wid, height=1, 
                             name=paste("vp.xlabels", x.name, sep=".")),
                 name=paste("text.xlabels", x.name, sep="."))
      allGrob <- gList(allGrob, 
                       gTree(children=gList(text.xlabels.sub),
                             vp=vp.xlabels,
                             name=paste("gTree.xlabels",
                                        x.name,
                                        sep=".")))
      
    }
  }    
  
  ## consider how to return all and could be redraw by grid.draw.
  grid.draw(allGrob)
  return(invisible(allGrob))
  }

# nathan_lawson new data 1(exclude no change ones)

# input.bed.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/BedWithLabel"

# input.bw.path <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/NR2F2kd_bigwig"
# dd <- 3000
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/NR2F2_exclude_no_change"
# null <- overLapWithOtherFeatures1(input.bed.dir,input.bw.path,output.file.dir,dd)

overLapWithOtherFeatures2 <- function(input.bed.dir,input.bw.path,output.file.dir,dd) {
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  file.name.4 <- list.files(input.bed.dir,pattern=".bed$",all.files = TRUE,full.names = TRUE,recursive = FALSE,include.dirs = TRUE)
  
  xx <- lapply(file.name.4, function(u){
    x <- readLines(u)
    x
  })
  
  yy <- lapply(xx, function(u){
    y <- as.data.frame(do.call(rbind,strsplit(u,split="\t")),stringsAsFactors=FALSE)
    y
  })
  
  yy[[2]] <- yy[[2]][-140,]
  yy[[4]] <- yy[[4]][-83,]
  
  names(yy) <- tools::file_path_sans_ext(basename(file.name.4))
  
  zz <- lapply(yy, function(u){
    if(dim(u)[2]<=3){
      z <- data.frame(u,paste0("no_",u[,1],"_",u[,2],"_",u[,3]))
    }
    else
    {
      u[grep("yes",u[,4]),4] <- paste0("yes_",u[grep("yes",u[,4]),1],"_",u[grep("yes",u[,4]),2],"_",u[grep("yes",u[,4]),3])
      
      u[-grep("yes",u[,4]),4] <- paste0("no_",u[-grep("yes",u[,4]),1],"_",u[-grep("yes",u[,4]),2],"_",u[-grep("yes",u[,4]),3])
      
      z <- u
    }
    
    colnames(z)=c("seqnames","start","end","label")
    
    z <- GRanges(z)
    z
  })
  
  generateMatrix4heatmap <- function(class_pattern,data,dd,input.bw.path,sortBy) {
    
    tss <- names(data)[grep(class_pattern,names(data),ignore.case = T)]
    print(tss)
    
    tss.bed <- data[which(names(data) %in% tss)]
    
    feature.center <- lapply(1:length(tss.bed), function(u,tss.bed,dd){
      x <- tss.bed[[u]]
      x
    },tss.bed,dd)
    
    names(feature.center) <- names(tss.bed)
    
    # get the regions of around feature middle +/- dd 
    feature.recentered.l <- lapply(1:length(tss.bed), function(u,tss.bed,dd){
      
      features <- tss.bed[[u]]
      
      wid <- width(features)
      feature.recentered <- feature.center <- features
      
      start(feature.center) <- start(features) + floor(wid/2)
      
      width(feature.center) <- 1
      
      start(feature.recentered) <- start(feature.center) - dd
      end(feature.recentered) <- end(feature.center) + dd
      
      feature.recentered
      
    },tss.bed,dd)
    names(feature.recentered.l) <- names(tss.bed)
    
    cvglists.l <- lapply(1:length(feature.recentered.l), function(u,feature.recentered.l,input.bw.path){
      
      files <- dir(input.bw.path, "bw")
      if(.Platform$OS.type != "windows"){
        cvglists <- sapply(file.path(input.bw.path,files), import, 
                           format="BigWig", 
                           which=feature.recentered.l[[u]], 
                           as="RleList")
      }else{## rtracklayer can not import bigWig files on Windows
        load(file.path(input.bw.path, "cvglist.rds"))
      }
      
    },feature.recentered.l,input.bw.path)
    
    tfType <- "TRC-NS"
    cvglists.ck <- lapply(1:length(cvglists.l), function(u,cvglists.l){
      z <- cvglists.l[[u]][grep(tfType,names(cvglists.l[[u]]))]
      z
    },cvglists.l)
    cvglists.ul.ck <- unlist(cvglists.ck)
    
    print(names(cvglists.ul.ck))
    
    n <- length(cvglists.ul.ck)
    
    sig1.4.heatmap.ck <- featureAlignedSignal2(cvglists.ck[[1]], feature.center[[1]], upstream=dd, downstream=dd)
    sig1.4.heatmap.log2.ck <- lapply(sig1.4.heatmap.ck, function(.ele) log2(.ele+1))
    names(sig1.4.heatmap.log2.ck) <- paste0(tfType,"_",seq(1,n))
    print(names(sig1.4.heatmap.log2.ck))
    
    tfType <- "shNR2F2"
    cvglists.tr <- lapply(1:length(cvglists.l), function(u,cvglists.l){
      z <- cvglists.l[[u]][grep(tfType,names(cvglists.l[[u]]))]
      z
    },cvglists.l)
    cvglists.ul.tr <- unlist(cvglists.tr)
    print(names(cvglists.ul.tr))
    n <- length(cvglists.ul.tr)
    
    sig1.4.heatmap.tr <- featureAlignedSignal2(cvglists.tr[[1]], feature.center[[1]], upstream=dd, downstream=dd)
    sig1.4.heatmap.log2.tr <- lapply(sig1.4.heatmap.tr, function(.ele) log2(.ele+1))
    names(sig1.4.heatmap.log2.tr) <- paste0(tfType,"_",seq(1,n))
    print(names(sig1.4.heatmap.log2.tr))
    
    dim(sig1.4.heatmap.log2.tr[[1]])
    
    #sortBy = "Trt"
    sig1.4.heatmap.log2.ck.tr.l <- lapply(1:length(sig1.4.heatmap.log2.tr),function(u,sig1.4.heatmap.log2.ck,sig1.4.heatmap.log2.tr,sortBy){
      if(sortBy=="Trt"){
        density.sum <- apply(sig1.4.heatmap.log2.tr[[u]],1,sum)
        sig1.4.heatmap.log2.tr.sorted<- sig1.4.heatmap.log2.tr[[u]][order(density.sum,decreasing = T),]        
        sig1.4.heatmap.log2.ck.sorted <- sig1.4.heatmap.log2.ck[[u]][match(rownames(sig1.4.heatmap.log2.tr.sorted),rownames(sig1.4.heatmap.log2.ck[[u]])),]
        sig1.4.heatmap.log2.ck.tr <- cbind(sig1.4.heatmap.log2.ck.sorted,sig1.4.heatmap.log2.tr.sorted)
      }else{
        density.sum <- apply(sig1.4.heatmap.log2.ck[[u]],1,sum)
        sig1.4.heatmap.log2.ck.sorted<- sig1.4.heatmap.log2.ck[[u]][order(density.sum,decreasing = T),]        
        sig1.4.heatmap.log2.tr.sorted <- sig1.4.heatmap.log2.tr[[u]][match(rownames(sig1.4.heatmap.log2.ck.sorted),rownames(sig1.4.heatmap.log2.tr[[u]])),]
        sig1.4.heatmap.log2.ck.tr <- cbind(sig1.4.heatmap.log2.ck.sorted,sig1.4.heatmap.log2.tr.sorted)
      }
      sig1.4.heatmap.log2.ck.tr
    },sig1.4.heatmap.log2.ck,sig1.4.heatmap.log2.tr,sortBy)
    
    sig1.4.heatmap.log2.ck.tr.l
    
  }
  
  heatmap4UpNoChangeDown <- function(class_patterns,file_name,data, dd, input.bw.path, output.file.dir,selected.genes=NULL,cols=colorpanel(30,"white","white","red")) {
    
    class_pattern <- class_patterns[1]
    sig1.4.heatmap.log2.ck.tr.Up <- generateMatrix4heatmap(class_pattern,data,dd,input.bw.path,sortBy="Trt")
    
    lbsUp <- as.character(mcols(data[[grep(class_pattern,names(data))]])$label)
    
    #print(names(sig1.4.heatmap.log2.ck.tr.Up))
    
    class_pattern <- class_patterns[2]
    sig1.4.heatmap.log2.ck.tr.Down <- generateMatrix4heatmap(class_pattern,data,dd,input.bw.path,sortBy="CK")
    
    lbsDown <- as.character(mcols(data[[grep(class_pattern,names(data))]])$label)
    #data[[grep(class_pattern,names(data))]]
    
    lbs <- c(lbsUp,lbsDown)
    
    if(is.null(selected.genes)) {selected.genes <- lbs[grep("yes",lbs)]}
    
    ht_list <- lapply(1:length(sig1.4.heatmap.log2.ck.tr.Up), function(u,sig1.4.heatmap.log2.ck.tr.Down,selected.genes){
      
      sig1.4.heatmap.log2.ck.tr.1 <- rbind(sig1.4.heatmap.log2.ck.tr.Up[[u]],sig1.4.heatmap.log2.ck.tr.Down[[u]])
      
      
      z.mat.0 <- t(scale(t(sig1.4.heatmap.log2.ck.tr.1), center=TRUE, scale=TRUE))
      
      print(rownames(z.mat.0))
      
      labels <- as.character(rownames(z.mat.0))
      selected.index <- match(selected.genes,labels)
      markers <- labels[selected.index]
      
      marker_idx <- selected.index
      
      print(labels[marker_idx])
      
      n.Up <- dim(sig1.4.heatmap.log2.ck.tr.Up[[u]])[1]
      n.Down <- dim(sig1.4.heatmap.log2.ck.tr.Down[[u]])[1]
      
      ht <- Heatmap(z.mat.0, name = paste0("z-score","-rep",u),
                    #col = rev(redblue(30))[-seq(35, 35)],  
                    #col = rev(redblue(30)),
                    # col = colorpanel(30,"white","white","red"),
                    col = cols,
                    show_row_name = FALSE,
                    cluster_columns = F,
                    cluster_rows = F,column_names_gp = gpar(fontsize = 6),column_split = rep(c("Control", "NR2F2KD"),c(100,100)),row_split = factor(rep(c("Up","Down"),c(n.Up,n.Down)),levels = c("Up","Down"))) + rowAnnotation(link = anno_mark(at = marker_idx, labels = markers),width = unit(5, "mm") +        max_text_width(markers))
      
      ht
      
    },sig1.4.heatmap.log2.ck.tr.Down,selected.genes)
    
    htList = ht_list[[1]] + ht_list[[2]] + ht_list[[3]]
    
    gb = grid.grabExpr(draw(htList))
    g <- plot_grid(gb)
    save_plot(file.path(output.file.dir,paste0(file_name,"_heatmap.png")),g,base_width = 25,base_height = 40)
  }
  
 
  
  class_patterns <- c("Art_K27acUp","Art_K27acDwn")
  file_name <- "Art_3_reps"
  data <- zz
  
  #[grep("Art",names(data))
  selected.genes <- c("yes_chr4_77370192_77371689")
  # for artery elements, let’s try red for the maximum density, ranging to white for the lowest;
  cols= colorpanel(100,"white","pink","red")
  
  # cols= colorpanel(100,"white","orange","red")
  
  heatmap4UpNoChangeDown(class_patterns,file_name,data, dd, input.bw.path, output.file.dir,cols=cols)
  
  # for vein let’s go with blue for max, white for minimum.
  
  class_patterns <- c("Vein_K27acUp","Vein_K27acDwn")
  file_name <- "Vein_3_reps"
  cols= colorpanel(100,"white","lightblue","blue")
  
  heatmap4UpNoChangeDown(class_patterns,file_name,data, dd, input.bw.path, output.file.dir,cols=cols)
  
}

fb <- function(n) {
  
  if(n == 0){x = 1}
  
  if(n ==1){x = 1}
  
  if(n >=2){ x = fb(n-1) + fb(n-2)}
  
  x  
}

featureAlignedSignal2 <- function(cvglists, feature.gr, 
                                  upstream, downstream, 
                                  n.tile=100, ...){
  stopifnot(is(feature.gr, "GRanges"))
  
  grWidr <- unique(width(feature.gr))
  if(missing(upstream) || missing(downstream)){
    if(length(grWidr)!=1){
      stop("The width of feature.gr is not identical.")
    }
  }else{
    if(!is.numeric(upstream) || !is.numeric(downstream)){
      stop("upstream and downstream must be integers")
    }
    if(upstream<0 || downstream<0){
      stop("upstream and downstream must be not less than 0")
    }
    upstream <- as.integer(upstream)
    downstream <- as.integer(downstream)
    if(length(grWidr)!=1 || any(grWidr!=1)){
      start(feature.gr) <- start(feature.gr)+floor(width(feature.gr)/2)
      width(feature.gr) <- 1
      warning("feature.gr is set to the center of feature.gr")
    }
    feature.gr <- promoters(feature.gr, upstream = upstream,
                            downstream = downstream + 1)
    grWidr <- unique(width(feature.gr))
  }
  if(any(start(feature.gr)<1)){
    warning("Some start position of the peaks are less than 1!",
            "They will be filtered.")
    feature.gr <- feature.gr[start(feature.gr)>0]
  }
  if(length(feature.gr)<2){
    stop("Length of feature.gr less than 2.")
  }
  if(inherits(cvglists, c("SimpleRleList", "RleList", "CompressedRleList"))){
    cvglistsName <- substitute(deparse(cvglists))
    cvglists <- list(cvglists)
    names(cvglists) <- cvglistsName
  }
  if(!is(cvglists, "list")){
    stop("cvglists must be a list of SimpleRleList or RleList")
  }
  cls <- sapply(cvglists, class)
  if(any(!cls %in% c("SimpleRleList", "RleList", "CompressedRleList")))
    stop("cvglists must be a list of SimpleRleList or RleList")
  seqLen <-  lapply(cvglists, function(.ele) 
    sapply(.ele, function(.e) sum(runLength(.e))))
  seqLen.keep <- table(unlist(sapply(seqLen, names)))==length(cvglists)
  seqLen <- seqLen[[1]][seqLen.keep]
  seqLen <- seqLen[!is.na(seqLen)]
  if(length(seqLen)>0){
    feature.gr.subset <- 
      feature.gr[seqnames(feature.gr) %in% names(seqLen)]
    if(any(end(feature.gr.subset)>seqLen[as.character(seqnames(feature.gr.subset))])){
      warning("Some end position of the peaks are out of bound!",
              "They will be filtered.")
      feature.gr <- 
        feature.gr.subset[end(feature.gr.subset) <=
                            seqLen[as.character(seqnames(feature.gr.subset))]]
    }
  }
  #feature.gr.bck <- feature.gr
  stopifnot(is.numeric(n.tile))
  n.tile <- round(n.tile)
  grL <- tile(feature.gr, n=n.tile)
  idx <- as.character(strand(feature.gr))=="-"
  if(sum(idx)>0) {
    grL.rev <- grL[idx]
    grL.rev.len <- lengths(grL.rev)
    grL.rev <- unlist(grL.rev, use.names = FALSE)
    grL.rev$oid <- rep(seq.int(length(grL[idx])), grL.rev.len)
    grL.rev <- rev(grL.rev)
    grL.rev.oid <- grL.rev$oid
    grL.rev$oid <- NULL
    grL.rev <- split(grL.rev, grL.rev.oid)
    grL[idx] <- as(grL.rev, "CompressedGRangesList")
    rm(grL.rev, grL.rev.len, grL.rev.oid)
  }
  grL.len <- lengths(grL)
  grL <- unlist(grL)
  grL$oid <- rep(seq_along(feature.gr), grL.len)
  grL$nid <- unlist(lapply(grL.len, seq.int))
  grL.len <- length(grL)
  grL.s <- split(grL, as.character(seqnames(grL)))
  grL <- unlist(grL.s) ## make sure grL and grL.s keep the same order
  seqn <- names(grL.s)
  seql <- seqlengths(feature.gr)
  seql <- seql[seqn]
  seql.f <- range(feature.gr)
  seql.f <- seql.f[match(seqn, seqnames(seql.f))]
  seql[is.na(seql)] <- width(seql.f)[is.na(seql)] + 1
  trimChar <- function(x, width){
    len <- nchar(x)
    if(len <= width) return(x)
    return(paste0(strtrim(x, width=width), "..."))
  }
  
  #rowname.feature.gr <- paste0(as.character(seqnames(feature.gr)), ":",
  #                             start(feature.gr), "-",
  #                             end(feature.gr))
  
  rowname.feature.gr <- as.character(mcols(feature.gr)$label)
  
  cov <- lapply(cvglists, function(.dat){
    .dat <- .dat[seqn[seqn %in% names(.dat)]]
    if(length(.dat)!=length(seqn)){
      warning(paste(seqn[!seqn %in% names(.dat)], collapse=", "), 
              ifelse(length(seqn[!seqn %in% names(.dat)])>1, " are", " is"), 
              " not in cvglists. seqlevels of cvglist are ", 
              trimChar(paste(names(.dat), collapse=", "), width=60))
      for(i in seqn[!seqn %in% names(.dat)]){
        .dat[[i]] <- Rle(0, seql[i])
      }
    }
    warn <- sapply(.dat, function(.ele){
      any(is.na(runValue(.ele)))
    })
    if(any(warn)){
      warning("cvglists contain NA values.")
    }
    .dat <- sapply(.dat, function(.ele){
      if(any(is.infinite(runValue(.ele)))){
        warning("cvglists contain infinite values. ", 
                "infinite value will be converted to NA.")
        runValue(.ele)[is.infinite(runValue(.ele))] <- NA
      }
      .ele
    })
    .dat <- .dat[seqn]
    vw <- Views(as(.dat, "RleList"), grL.s)
    vm <- viewMeans(vw)
    vm <- unlist(vm)
    stopifnot(length(vm)==grL.len)
    mm <- matrix(0, nrow=length(feature.gr), ncol=n.tile)
    mm[grL$oid+length(feature.gr)*(grL$nid-1)] <- vm
    rownames(mm) <- rowname.feature.gr
    mm
  })
  
  cov
}

ivanPlay <- function(){
  cat("hi my name is Ivan Yan!")
}

# input.bed.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/2019_1_7_forHeatMapsAndVennDiagrams"
# input.bw.path <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/ERGChipSeq_bigwig/bigwig"
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/Results_4_VeNR2F2repsFinal"
# dd <- 3000

overLapWithOtherFeatures3 <- function(input.bed.dir,input.bw.path,output.file.dir,dd) {
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  file.name.4 <- list.files(input.bed.dir,pattern=".bed$",all.files = TRUE,full.names = TRUE,recursive = FALSE,include.dirs = TRUE)
  
  xx <- lapply(file.name.4, function(u){
    x <- readLines(u)
    x
  })
  
  yy <- lapply(xx, function(u){
    y <- as.data.frame(do.call(rbind,strsplit(u,split="\t")),stringsAsFactors=FALSE)
    y
  })
  
  yy[[2]] <- yy[[2]][-140,]
  yy[[4]] <- yy[[4]][-83,]
  
  names(yy) <- tools::file_path_sans_ext(basename(file.name.4))
  
  zz <- lapply(yy, function(u){
    z <- u[,c(1:3)]
    
    z <- data.frame(u[,c(1:3)],paste0(as.character(u[,1]), ":",u[,2], "-",u[,3]))
    
    colnames(z)=c("seqnames","start","end","label")
    
    z <- GRanges(z)
    z
  })
  
  generateMatrix4heatmap <- function(class_pattern,data,dd,input.bw.path,sortBy) {
    
    # class_pattern <- "VeNR2F2repsFinal"
    #  file_name <- "Art_3_reps"
    #  data <- zz
    
    tss <- names(data)[grep(class_pattern,names(data),ignore.case = T)]
    print(tss)
    
    tss.bed <- data[which(names(data) %in% tss)]
    
    feature.center <- lapply(1:length(tss.bed), function(u,tss.bed,dd){
      x <- tss.bed[[u]]
      x
    },tss.bed,dd)
    
    names(feature.center) <- names(tss.bed)
    
    # get the regions of around feature middle +/- dd 
    feature.recentered.l <- lapply(1:length(tss.bed), function(u,tss.bed,dd){
      
      features <- tss.bed[[u]]
      
      wid <- width(features)
      feature.recentered <- feature.center <- features
      
      start(feature.center) <- start(features) + floor(wid/2)
      
      width(feature.center) <- 1
      
      start(feature.recentered) <- start(feature.center) - dd
      end(feature.recentered) <- end(feature.center) + dd
      
      feature.recentered
      
    },tss.bed,dd)
    names(feature.recentered.l) <- names(tss.bed)
    
    cvglists.l <- lapply(1:length(feature.recentered.l), function(u,feature.recentered.l,input.bw.path){
      
      files <- dir(input.bw.path, "bw")
      if(.Platform$OS.type != "windows"){
        cvglists <- sapply(file.path(input.bw.path,files), import, 
                           format="BigWig", 
                           which=feature.recentered.l[[u]], 
                           as="RleList")
      }else{## rtracklayer can not import bigWig files on Windows
        load(file.path(input.bw.path, "cvglist.rds"))
      }
      
    },feature.recentered.l,input.bw.path)
    
    tfType <- "V42ERG_IN"
    cvglists.ck <- lapply(1:length(cvglists.l), function(u,cvglists.l){
      z <- cvglists.l[[u]][grep(tfType,names(cvglists.l[[u]]))]
      z
    },cvglists.l)
    cvglists.ul.ck <- unlist(cvglists.ck)
    
    print(names(cvglists.ul.ck))
    
    n <- length(cvglists.ul.ck)
    
    sig1.4.heatmap.ck <- featureAlignedSignal2(cvglists.ck[[1]], feature.center[[1]], upstream=dd, downstream=dd)
    sig1.4.heatmap.log2.ck <- lapply(sig1.4.heatmap.ck, function(.ele) log2(.ele+1))
    names(sig1.4.heatmap.log2.ck) <- paste0(tfType,seq(1,n))
    print(names(sig1.4.heatmap.log2.ck))
    
    tfType <- "V42ERG_IP"
    cvglists.tr <- lapply(1:length(cvglists.l), function(u,cvglists.l){
      z <- cvglists.l[[u]][grep(tfType,names(cvglists.l[[u]]))]
      z
    },cvglists.l)
    cvglists.ul.tr <- unlist(cvglists.tr)
    print(names(cvglists.ul.tr))
    n <- length(cvglists.ul.tr)
    
    sig1.4.heatmap.tr <- featureAlignedSignal2(cvglists.tr[[1]], feature.center[[1]], upstream=dd, downstream=dd)
    sig1.4.heatmap.log2.tr <- lapply(sig1.4.heatmap.tr, function(.ele) log2(.ele+1))
    names(sig1.4.heatmap.log2.tr) <- paste0(tfType,seq(1,n))
    print(names(sig1.4.heatmap.log2.tr))
    
    dim(sig1.4.heatmap.log2.tr[[1]])
    
    sortBy = "Trt"
    sig1.4.heatmap.log2.ck.tr.l <- lapply(1:length(sig1.4.heatmap.log2.tr),function(u,sig1.4.heatmap.log2.ck,sig1.4.heatmap.log2.tr,sortBy){
      if(sortBy=="Trt"){
        density.sum <- apply(sig1.4.heatmap.log2.tr[[u]],1,sum)
        sig1.4.heatmap.log2.tr.sorted<- sig1.4.heatmap.log2.tr[[u]][order(density.sum,decreasing = T),]        
        sig1.4.heatmap.log2.ck.sorted <- sig1.4.heatmap.log2.ck[[u]][match(rownames(sig1.4.heatmap.log2.tr.sorted),rownames(sig1.4.heatmap.log2.ck[[u]])),]
        sig1.4.heatmap.log2.ck.tr <- cbind(sig1.4.heatmap.log2.ck.sorted,sig1.4.heatmap.log2.tr.sorted)
      }else{
        density.sum <- apply(sig1.4.heatmap.log2.ck[[u]],1,sum)
        sig1.4.heatmap.log2.ck.sorted<- sig1.4.heatmap.log2.ck[[u]][order(density.sum,decreasing = T),]        
        sig1.4.heatmap.log2.tr.sorted <- sig1.4.heatmap.log2.tr[[u]][match(rownames(sig1.4.heatmap.log2.ck.sorted),rownames(sig1.4.heatmap.log2.tr[[u]])),]
        sig1.4.heatmap.log2.ck.tr <- cbind(sig1.4.heatmap.log2.ck.sorted,sig1.4.heatmap.log2.tr.sorted)
      }
      sig1.4.heatmap.log2.ck.tr
    },sig1.4.heatmap.log2.ck,sig1.4.heatmap.log2.tr,sortBy)
    
    sig1.4.heatmap.log2.ck.tr.l
    
  }
  
  heatmap4UpNoChangeDown <- function(class_patterns,file_name,data, dd, input.bw.path, output.file.dir,selected.genes=NULL,cols=colorpanel(30,"white","white","red")) {
    
    class_patterns <- c("VeNR2F2repsFinal")
    
    class_patterns <- c("VeERGrepsFINAL")
    
    class_pattern <- class_patterns[1]
    sig1.4.heatmap.log2.ck.tr.Up <- generateMatrix4heatmap(class_pattern,data,dd,input.bw.path,sortBy="Trt")
    
    lbsUp <- as.character(mcols(data[[grep(class_pattern,names(data))]])$label)
    
    #print(names(sig1.4.heatmap.log2.ck.tr.Up))
    
    class_pattern <- class_patterns[2]
    sig1.4.heatmap.log2.ck.tr.Down <- generateMatrix4heatmap(class_pattern,data,dd,input.bw.path,sortBy="CK")
    
    lbsDown <- as.character(mcols(data[[grep(class_pattern,names(data))]])$label)
    #data[[grep(class_pattern,names(data))]]
    
    lbs <- c(lbsUp,lbsDown)
    
    if(is.null(selected.genes)) {selected.genes <- lbs[grep("yes",lbs)]}
    
    cols <- colorpanel(100,"white","lightblue","blue")
    
    ht_list <- lapply(1:length(sig1.4.heatmap.log2.ck.tr.Up), function(u,cols){
      
      #u <- 1
      sig1.4.heatmap.log2.ck.tr.1 <- sig1.4.heatmap.log2.ck.tr.Up[[u]]
      
      
      z.mat.0 <- t(scale(t(sig1.4.heatmap.log2.ck.tr.1), center=TRUE, scale=TRUE))
      
      print(rownames(z.mat.0))
      
      ht <- Heatmap(z.mat.0, name = paste0("z-score","-rep",u),
                    #col = rev(redblue(30))[-seq(35, 35)],  
                    #col = rev(redblue(30)),
                    # col = colorpanel(30,"white","white","red"),
                    col = cols,
                    show_row_name = FALSE,
                    cluster_columns = F,
                    cluster_rows = F,column_names_gp = gpar(fontsize = 6),column_split = rep(c(paste0("V42ERG_IN",u), paste0("V42ERG_IP",u)),c(100,100)))
      
      ht
      
    },cols)
    
    htList = ht_list[[1]] + ht_list[[2]]
    
    gb = grid.grabExpr(draw(htList))
    g <- plot_grid(gb)
    save_plot(file.path(output.file.dir,paste0(file_name,"_heatmap.png")),g,base_width = 25,base_height = 40)
  }
  
  library(ComplexHeatmap)
  library(gplots)
  library(cowplot)
  
  class_patterns <- c("Art_K27acUp","Art_K27acDwn")
  file_name <- "Art_3_reps"
  data <- zz
  
  #[grep("Art",names(data))
  selected.genes <- c("yes_chr4_77370192_77371689")
  # for artery elements, let’s try red for the maximum density, ranging to white for the lowest;
  cols= colorpanel(100,"white","pink","red")
  
  # cols= colorpanel(100,"white","orange","red")
  
  heatmap4UpNoChangeDown(class_patterns,file_name,data, dd, input.bw.path, output.file.dir,cols=cols)
  
  # for vein let’s go with blue for max, white for minimum.
  
  class_patterns <- c("Vein_K27acUp","Vein_K27acDwn")
  file_name <- "Vein_3_reps"
  cols= colorpanel(100,"white","lightblue","blue")
  
  heatmap4UpNoChangeDown(class_patterns,file_name,data, dd, input.bw.path, output.file.dir,cols=cols)
  
  
  output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/Results_4_VeNR2F2repsFinal/VennByRegions"
  
  peak.index <- c(3,5)
  name <- c("VeERGrepsFINAL","VeNR2F2repsFinal")
  getCount4Venn(zz,peak.index,name,output.file.dir)
  
  grl <- zz[peak.index]
  pdf(file.path(output.file.dir,paste0(paste(name,collapse ="_"),".pdf")))
  makeVennDiagram(grl, NameOfPeaks=name,fill = c("blue", "green"),cat.col = c("blue", "green"),euler.d = TRUE,cex = 1,cat.dist=0.01,totalTest=56967)
  dev.off()
  
  getCount4Venn2(zz,peak.index,name,output.file.dir)
  
  # VeERGrepsFINAL.bed
  # VeNR2F2repsFinal.bed
  # ArteryAll.bed
  
  peak.index <- c(3,5,1)
  name <- c("VeERGrepsFINAL","VeNR2F2repsFinal","ArteryAll")
  getCount4Venn(zz,peak.index,name,output.file.dir)
  
  grl <- zz[peak.index]
  pdf(file.path(output.file.dir,paste0(paste(name,collapse ="_"),".pdf")))
  makeVennDiagram(grl, NameOfPeaks=name,fill = c("blue", "green","red"),cat.col = c("blue", "green","red"),euler.d = TRUE,cex = 1,cat.dist=0.01,totalTest=56967)
  dev.off() 
  
  # VeERGrepsFINAL.bed
  # VeNR2F2repsFinal.bed
  # CommonAll.bed				
  peak.index <- c(3,5,2)
  name <- c("VeERGrepsFINAL","VeNR2F2repsFinal","CommonAll")
  getCount4Venn(zz,peak.index,name,output.file.dir)
  
  
  grl <- zz[peak.index]
  pdf(file.path(output.file.dir,paste0(paste(name,collapse ="_"),".pdf")))
  makeVennDiagram(grl, NameOfPeaks=name,fill = c("blue", "green","red"),cat.col = c("blue", "green","red"),euler.d = TRUE,cex = 1,cat.dist=0.01,totalTest=56967)
  dev.off() 
  
  
  # VeERGrepsFINAL.bed
  # VeNR2F2repsFinal.bed
  # VeinAll.bed
  peak.index <- c(3,5,4)
  name <- c("VeERGrepsFINAL","VeNR2F2repsFinal","VeinAll")
  getCount4Venn(zz,peak.index,name,output.file.dir)
  
  grl <- zz[peak.index]
  pdf(file.path(output.file.dir,paste0(paste(name,collapse ="_"),".pdf")))
  makeVennDiagram(grl, NameOfPeaks=name,fill = c("blue", "green","red"),cat.col = c("blue", "green","red"),euler.d = TRUE,cex = 1,cat.dist=0.01,totalTest=56967)
  dev.off() 
  
}

# input.bed.dir <- "~/Aimin/project/umw_nathan_lawson/deep_seq_data/CutAndRun/BedFilesForPlotsAndHeatMaps"

# input.bw.path <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/OutPutCutAndRun"

# input.bw.path <- "~/Aimin/nl/umw_nathan_lawson/NR2F2bigwigs"

# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/Results_4_CutAndRun"

# dd <- 10000

overLapWithOtherFeatures4 <- function(input.bed.dir,input.bw.path,output.file.dir,dd) {
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  file.name.4 <- list.files(input.bed.dir,pattern=".bed$",all.files = TRUE,full.names = TRUE,recursive = FALSE,include.dirs = TRUE)
  
  xx <- lapply(file.name.4, function(u){
    x <- readLines(u)
    x
  })
  
  yy <- lapply(xx, function(u){
    y <- as.data.frame(do.call(rbind,strsplit(u,split="\t")),stringsAsFactors=FALSE)
    y
  })
  
  yy[[4]] <- yy[[4]][-1,]
  yy[[7]] <- yy[[7]][-1,]
  yy[[10]] <- yy[[10]][-1,]

  yy[[3]] <- yy[[3]][-1,]
  yy[[6]] <- yy[[6]][-1,]
  yy[[9]] <- yy[[9]][-1,]
  
  # ArteryAll.bed, CommonAll.bed,VeinAll.bed
  
  # ArteryTSS.bed,CommonTSS.bed, VeinTSS.bed
  
  # ArteryEnh.bed,CommonEnh.bed, VeinEnh.bed
  
  names(yy) <- tools::file_path_sans_ext(basename(file.name.4))
  
  zz <- lapply(yy[c(2,5,8,4,7,10,3,6,9)], function(u){
    #z <- u[,c(1:3)]
    
    z <- data.frame(u[,c(1:3)],paste0(as.character(u[,1]), ":",u[,2], "-",u[,3]))
    
    colnames(z)=c("seqnames","start","end","label")
    
    z <- GRanges(z)
    z
  })
  
  generateMatrix4heatmap <- function(class_pattern,data,dd,input.bw.path,sortBy) {
    
    # class_pattern <- "ArteryAll"
    # file_name <- "HDAC1"
    # data <- zz
    
    tss <- names(data)[grep(class_pattern,names(data),ignore.case = T)]
    print(tss)
    
    tss.bed <- data[which(names(data) %in% tss)]
    
    feature.center <- lapply(1:length(tss.bed), function(u,tss.bed,dd){
      x <- tss.bed[[u]]
      x
    },tss.bed,dd)
    
    names(feature.center) <- names(tss.bed)
    
    # get the regions of around feature middle +/- dd 
    feature.recentered.l <- lapply(1:length(tss.bed), function(u,tss.bed,dd){
      
      features <- tss.bed[[u]]
      
      wid <- width(features)
      feature.recentered <- feature.center <- features
      
      start(feature.center) <- start(features) + floor(wid/2)
      
      width(feature.center) <- 1
      
      start(feature.recentered) <- start(feature.center) - dd
      end(feature.recentered) <- end(feature.center) + dd
      
      feature.recentered
      
    },tss.bed,dd)
    names(feature.recentered.l) <- names(tss.bed)
    
    cvglists.l <- lapply(1:length(feature.recentered.l), function(u,feature.recentered.l,input.bw.path){
      
      files <- dir(input.bw.path, "bw")
      if(.Platform$OS.type != "windows"){
        cvglists <- sapply(file.path(input.bw.path,files), import, 
                           format="BigWig", 
                           which=feature.recentered.l[[u]], 
                           as="RleList")
      }else{## rtracklayer can not import bigWig files on Windows
        load(file.path(input.bw.path, "cvglist.rds"))
      }
      
    },feature.recentered.l,input.bw.path)
    
    tfType <- "HDAC1"
    cvglists.ck <- lapply(1:length(cvglists.l), function(u,cvglists.l){
      z <- cvglists.l[[u]][grep(tfType,names(cvglists.l[[u]]))]
      z
    },cvglists.l)
    cvglists.ul.ck <- unlist(cvglists.ck)
    
    print(names(cvglists.ul.ck))
    
    n <- length(cvglists.ul.ck)
    
    sig1.4.heatmap.ck <- featureAlignedSignal2(cvglists.ck[[1]], feature.center[[1]], upstream=dd, downstream=dd)
    sig1.4.heatmap.log2.ck <- lapply(sig1.4.heatmap.ck, function(.ele) log2(.ele+1))
    names(sig1.4.heatmap.log2.ck) <- paste0(tfType,seq(1,n))
    print(names(sig1.4.heatmap.log2.ck))
    
    out.4 <- featureAlignedDistribution(cvglists.ck[[1]], 
                                        feature.center[[1]],upstream=dd, downstream=dd,
                                        zeroAt=.5,type="l", 
                                        ylab="Averaged coverage",add=TRUE)
    
    x.pos <- seq(-dd,dd,length.out = dim(out)[1])
    out.m <- data.frame(x.pos=x.pos,out.4)
    
    write.table(out.m,file = file.path(output.file.dir,paste0(class_pattern,"_Density_matrix_short",".txt")),
                append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T)
    
    
    tfType <- "V42ERG_IP"
    cvglists.tr <- lapply(1:length(cvglists.l), function(u,cvglists.l){
      z <- cvglists.l[[u]][grep(tfType,names(cvglists.l[[u]]))]
      z
    },cvglists.l)
    cvglists.ul.tr <- unlist(cvglists.tr)
    print(names(cvglists.ul.tr))
    n <- length(cvglists.ul.tr)
    
    sig1.4.heatmap.tr <- featureAlignedSignal2(cvglists.tr[[1]], feature.center[[1]], upstream=dd, downstream=dd)
    sig1.4.heatmap.log2.tr <- lapply(sig1.4.heatmap.tr, function(.ele) log2(.ele+1))
    names(sig1.4.heatmap.log2.tr) <- paste0(tfType,seq(1,n))
    print(names(sig1.4.heatmap.log2.tr))
    
    dim(sig1.4.heatmap.log2.tr[[1]])
    
    sortBy = "Trt"
    sig1.4.heatmap.log2.ck.tr.l <- lapply(1:length(sig1.4.heatmap.log2.tr),function(u,sig1.4.heatmap.log2.ck,sig1.4.heatmap.log2.tr,sortBy){
      if(sortBy=="Trt"){
        density.sum <- apply(sig1.4.heatmap.log2.tr[[u]],1,sum)
        sig1.4.heatmap.log2.tr.sorted<- sig1.4.heatmap.log2.tr[[u]][order(density.sum,decreasing = T),]        
        sig1.4.heatmap.log2.ck.sorted <- sig1.4.heatmap.log2.ck[[u]][match(rownames(sig1.4.heatmap.log2.tr.sorted),rownames(sig1.4.heatmap.log2.ck[[u]])),]
        sig1.4.heatmap.log2.ck.tr <- cbind(sig1.4.heatmap.log2.ck.sorted,sig1.4.heatmap.log2.tr.sorted)
      }else{
        density.sum <- apply(sig1.4.heatmap.log2.ck[[u]],1,sum)
        sig1.4.heatmap.log2.ck.sorted<- sig1.4.heatmap.log2.ck[[u]][order(density.sum,decreasing = T),]        
        sig1.4.heatmap.log2.tr.sorted <- sig1.4.heatmap.log2.tr[[u]][match(rownames(sig1.4.heatmap.log2.ck.sorted),rownames(sig1.4.heatmap.log2.tr[[u]])),]
        sig1.4.heatmap.log2.ck.tr <- cbind(sig1.4.heatmap.log2.ck.sorted,sig1.4.heatmap.log2.tr.sorted)
      }
      sig1.4.heatmap.log2.ck.tr
    },sig1.4.heatmap.log2.ck,sig1.4.heatmap.log2.tr,sortBy)
    
    sig1.4.heatmap.log2.ck.tr.l
    
  }
  
  heatmap4UpNoChangeDown <- function(class_patterns,file_name,data, dd, input.bw.path, output.file.dir,selected.genes=NULL,cols=colorpanel(30,"white","white","red")) {
    
    class_patterns <- c("VeNR2F2repsFinal")
    
    class_patterns <- c("VeERGrepsFINAL")
    
    class_pattern <- class_patterns[1]
    sig1.4.heatmap.log2.ck.tr.Up <- generateMatrix4heatmap(class_pattern,data,dd,input.bw.path,sortBy="Trt")
    
    lbsUp <- as.character(mcols(data[[grep(class_pattern,names(data))]])$label)
    
    #print(names(sig1.4.heatmap.log2.ck.tr.Up))
    
    class_pattern <- class_patterns[2]
    sig1.4.heatmap.log2.ck.tr.Down <- generateMatrix4heatmap(class_pattern,data,dd,input.bw.path,sortBy="CK")
    
    lbsDown <- as.character(mcols(data[[grep(class_pattern,names(data))]])$label)
    #data[[grep(class_pattern,names(data))]]
    
    lbs <- c(lbsUp,lbsDown)
    
    if(is.null(selected.genes)) {selected.genes <- lbs[grep("yes",lbs)]}
    
    cols <- colorpanel(100,"white","lightblue","blue")
    
    ht_list <- lapply(1:length(sig1.4.heatmap.log2.ck.tr.Up), function(u,cols){
      
      #u <- 1
      sig1.4.heatmap.log2.ck.tr.1 <- sig1.4.heatmap.log2.ck.tr.Up[[u]]
      
      
      z.mat.0 <- t(scale(t(sig1.4.heatmap.log2.ck.tr.1), center=TRUE, scale=TRUE))
      
      print(rownames(z.mat.0))
      
      ht <- Heatmap(z.mat.0, name = paste0("z-score","-rep",u),
                    #col = rev(redblue(30))[-seq(35, 35)],  
                    #col = rev(redblue(30)),
                    # col = colorpanel(30,"white","white","red"),
                    col = cols,
                    show_row_name = FALSE,
                    cluster_columns = F,
                    cluster_rows = F,column_names_gp = gpar(fontsize = 6),column_split = rep(c(paste0("V42ERG_IN",u), paste0("V42ERG_IP",u)),c(100,100)))
      
      ht
      
    },cols)
    
    htList = ht_list[[1]] + ht_list[[2]]
    
    gb = grid.grabExpr(draw(htList))
    g <- plot_grid(gb)
    save_plot(file.path(output.file.dir,paste0(file_name,"_heatmap.png")),g,base_width = 25,base_height = 40)
  }
  
  library(ComplexHeatmap)
  library(gplots)
  library(cowplot)
  
  class_patterns <- c("Art_K27acUp","Art_K27acDwn")
  file_name <- "Art_3_reps"
  data <- zz
  
  #[grep("Art",names(data))
  selected.genes <- c("yes_chr4_77370192_77371689")
  # for artery elements, let’s try red for the maximum density, ranging to white for the lowest;
  cols= colorpanel(100,"white","pink","red")
  
  # cols= colorpanel(100,"white","orange","red")
  
  heatmap4UpNoChangeDown(class_patterns,file_name,data, dd, input.bw.path, output.file.dir,cols=cols)
  
  # for vein let’s go with blue for max, white for minimum.
  
  class_patterns <- c("Vein_K27acUp","Vein_K27acDwn")
  file_name <- "Vein_3_reps"
  cols= colorpanel(100,"white","lightblue","blue")
  
  heatmap4UpNoChangeDown(class_patterns,file_name,data, dd, input.bw.path, output.file.dir,cols=cols)
  
  
  output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/Results_4_VeNR2F2repsFinal/VennByRegions"
  
  peak.index <- c(3,5)
  name <- c("VeERGrepsFINAL","VeNR2F2repsFinal")
  getCount4Venn(zz,peak.index,name,output.file.dir)
  
  grl <- zz[peak.index]
  pdf(file.path(output.file.dir,paste0(paste(name,collapse ="_"),".pdf")))
  makeVennDiagram(grl, NameOfPeaks=name,fill = c("blue", "green"),cat.col = c("blue", "green"),euler.d = TRUE,cex = 1,cat.dist=0.01,totalTest=56967)
  dev.off()
  
  getCount4Venn2(zz,peak.index,name,output.file.dir)
  
  # VeERGrepsFINAL.bed
  # VeNR2F2repsFinal.bed
  # ArteryAll.bed
  
  peak.index <- c(3,5,1)
  name <- c("VeERGrepsFINAL","VeNR2F2repsFinal","ArteryAll")
  getCount4Venn(zz,peak.index,name,output.file.dir)
  
  grl <- zz[peak.index]
  pdf(file.path(output.file.dir,paste0(paste(name,collapse ="_"),".pdf")))
  makeVennDiagram(grl, NameOfPeaks=name,fill = c("blue", "green","red"),cat.col = c("blue", "green","red"),euler.d = TRUE,cex = 1,cat.dist=0.01,totalTest=56967)
  dev.off() 
  
  # VeERGrepsFINAL.bed
  # VeNR2F2repsFinal.bed
  # CommonAll.bed				
  peak.index <- c(3,5,2)
  name <- c("VeERGrepsFINAL","VeNR2F2repsFinal","CommonAll")
  getCount4Venn(zz,peak.index,name,output.file.dir)
  
  
  grl <- zz[peak.index]
  pdf(file.path(output.file.dir,paste0(paste(name,collapse ="_"),".pdf")))
  makeVennDiagram(grl, NameOfPeaks=name,fill = c("blue", "green","red"),cat.col = c("blue", "green","red"),euler.d = TRUE,cex = 1,cat.dist=0.01,totalTest=56967)
  dev.off() 
  
  # VeERGrepsFINAL.bed
  # VeNR2F2repsFinal.bed
  # VeinAll.bed
  peak.index <- c(3,5,4)
  name <- c("VeERGrepsFINAL","VeNR2F2repsFinal","VeinAll")
  getCount4Venn(zz,peak.index,name,output.file.dir)
  
  grl <- zz[peak.index]
  pdf(file.path(output.file.dir,paste0(paste(name,collapse ="_"),".pdf")))
  makeVennDiagram(grl, NameOfPeaks=name,fill = c("blue", "green","red"),cat.col = c("blue", "green","red"),euler.d = TRUE,cex = 1,cat.dist=0.01,totalTest=56967)
  dev.off() 
  
}

overLapWithOtherFeatures5 <- function(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType) {
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  file.name.4 <- list.files(input.bed.dir,pattern=".bed$",all.files = TRUE,full.names = TRUE,recursive = FALSE,include.dirs = TRUE)
  
  xx <- lapply(file.name.4, function(u){
    x <- readLines(u)
    x
  })
  
  yy <- lapply(xx, function(u){
    y <- as.data.frame(do.call(rbind,strsplit(u,split="\t")),stringsAsFactors=FALSE)
    y
  })
  
  yy[[4]] <- yy[[4]][-1,]
  yy[[7]] <- yy[[7]][-1,]
  yy[[10]] <- yy[[10]][-1,]
  
  yy[[3]] <- yy[[3]][-1,]
  yy[[6]] <- yy[[6]][-1,]
  yy[[9]] <- yy[[9]][-1,]
  
  # ArteryAll.bed, CommonAll.bed,VeinAll.bed
  
  # ArteryTSS.bed,CommonTSS.bed, VeinTSS.bed
  
  # ArteryEnh.bed,CommonEnh.bed, VeinEnh.bed
  
  names(yy) <- tools::file_path_sans_ext(basename(file.name.4))
  
  zz <- lapply(yy[c(2,5,8,4,7,10,3,6,9)], function(u){
    #z <- u[,c(1:3)]
    
    z <- data.frame(u[,c(1:3)],paste0(as.character(u[,1]), ":",u[,2], "-",u[,3]))
    
    colnames(z)=c("seqnames","start","end","label")
    
    z <- GRanges(z)
    z
  })
  
  data <- zz
  
  class_pattern <- "ArteryAll"
  
  tss <- names(data)[grep(class_pattern,names(data),ignore.case = T)]
  print(tss)
  
  tss.bed <- data[which(names(data) %in% tss)]
  
  feature.center <- lapply(1:length(tss.bed), function(u,tss.bed,dd){
    x <- tss.bed[[u]]
    x
  },tss.bed,dd)
  
  names(feature.center) <- names(tss.bed)
  
  # get the regions of around feature middle +/- dd 
  feature.recentered.l <- lapply(1:length(tss.bed), function(u,tss.bed,dd){
    
    features <- tss.bed[[u]]
    
    wid <- width(features)
    feature.recentered <- feature.center <- features
    
    start(feature.center) <- start(features) + floor(wid/2)
    
    width(feature.center) <- 1
    
    start(feature.recentered) <- start(feature.center) - dd
    end(feature.recentered) <- end(feature.center) + dd
    
    feature.recentered
    
  },tss.bed,dd)
  
  names(feature.recentered.l) <- names(tss.bed)

  input.bw.path <- "~/Aimin/nl/umw_nathan_lawson/NR2F2bigwigs"
  
  cvglists.l <- lapply(1:length(feature.recentered.l), function(u,feature.recentered.l,input.bw.path){
    
    files <- dir(input.bw.path, "bigwig")
    if(.Platform$OS.type != "windows"){
      cvglists <- sapply(file.path(input.bw.path,files), import, 
                         format="BigWig", 
                         which=feature.recentered.l[[u]], 
                         as="RleList")
    }else{## rtracklayer can not import bigWig files on Windows
      load(file.path(input.bw.path, "cvglist.rds"))
    }
    
  },feature.recentered.l,input.bw.path)
  
  #cvglists.l <- HDAC1_10000_ArteryTSS_CommonTSS_VeinTSS_cvglist
  tfType <- "V29NR2F2_"
  
  grep("V29NR2F2_",names(unlist(cvglists.l)))
  
  print(length(cvglists.l))
  print(length(tss.bed))
  
  #x_name_0 <- tools::file_path_sans_ext(basename(names(cvglists.l[[1]])))  
  #x_name_start <- str_locate_all(x_name_0,"_")[[1]][2,2] + 1
  #x_name_end <- str_locate_all(x_name_0,"_")[[1]][3,2] - 1
  #x_name <- paste0(unique(str_sub(x_name_0,x_name_start,x_name_end)),"_",dd,"_",paste(tss,collapse = "_"))
  
  x_name <- paste0(tfType,"_",dd,"_",paste(tss,collapse = "_"))
  
  saveRDS(cvglists.l,file = file.path(output.file.dir,paste0(x_name,"_cvglist.rds")))
  
  
  
  #lapply(1:length(tss.bed),function(u){
  
  #Class1 <- (cvglists.ul[[1]]+cvglists.ul[[2]])/2
  
  #Class2 <- (cvglists.ul[[3]]+cvglists.ul[[4]])/2
  
  #Class3 <- (cvglists.ul[[5]]+cvglists.ul[[6]])/2
  
  #Endo <- (cvglists.ul[[7]]+cvglists.ul[[8]])/2
  
  #}
  
  #cvglists.ave <- list(Class1=Class1,Class2=Class2,Class3=Class3, Endo=Endo)
  
  cvglists.l <- lapply(1:length(cvglists.l), function(u,cvglists.l){
    
    z <- cvglists.l[[u]][grep(tfType,names(cvglists.l[[u]]))]
    z
  },cvglists.l)
  
  #names(cvglists.l[[4]])
  
  cvglists.ul <- unlist(cvglists.l)
  
  #cvglists.ul <- cvglists.ul[grep(tfType,names(cvglists.ul))]
  
  #grep(tfType,names(unlist(cvglists.l)))
  
  n <- length(cvglists.ul)/length(tss.bed)
  
  if(n==2){
    cvglists.ave <- lapply(1:length(tss.bed),function(u,cvglists.ul){
      
      Class1 <- (cvglists.ul[[2*u-1]]+cvglists.ul[[2*u]])/n
      
      Class1
      
    },cvglists.ul)
  }
  
  if(n==3){
    cvglists.ave <- lapply(1:length(tss.bed),function(u,cvglists.ul){
      
      Class1 <- (cvglists.ul[[3*u-2]]+cvglists.ul[[3*u-1]]+cvglists.ul[[3*u]])/n
      
      Class1
      
    },cvglists.ul)
  }
  
  #cvglists.ave <- list(Class1=Class1,Class2=Class2,Class3=Class3, Endo=Endo)
  
  names(cvglists.ave) <- names(tss.bed)
  
  
  
  
  # sig1 <- featureAlignedSignal(cvglists.ave, feature.center[[1]], upstream=dd, downstream=dd) 
  # sigs.log2 <- lapply(sig1, function(.ele) log2(.ele+1))   
  # out.1 <- featureAlignedDistribution(sigs.log2, 
  #                                   feature.center[[1]],upstream=dd, downstream=dd,
  #                                   zeroAt=.5,type="l", 
  #                                   ylab="Averaged coverage")
  # 
  # 
  # sig1 <- featureAlignedSignal(cvglists.ave, feature.center[[2]], upstream=dd, downstream=dd) 
  # sigs.log2 <- lapply(sig1, function(.ele) log2(.ele+1))
  # out.2 <- featureAlignedDistribution(sigs.log2, 
  #                                   feature.center[[2]],upstream=dd, downstream=dd,
  #                                   zeroAt=.5,type="l", 
  #                                   ylab="Averaged coverage",add=TRUE)
  # 
  # sig1 <- featureAlignedSignal(cvglists.ave, feature.center[[3]], upstream=dd, downstream=dd) 
  # sigs.log2 <- lapply(sig1, function(.ele) log2(.ele+1))
  # out.3 <- featureAlignedDistribution(sigs.log2, 
  #                                   feature.center[[3]],upstream=dd, downstream=dd,
  #                                   zeroAt=.5,type="l", 
  #                                   ylab="Averaged coverage",add=TRUE)
  # 
  # 
  # sig1 <- featureAlignedSignal(cvglists.ave, feature.center[[4]], upstream=dd, downstream=dd) 
  # sigs.log2 <- lapply(sig1, function(.ele) log2(.ele+1))
  # out.4 <- featureAlignedDistribution(sigs.log2, 
  #                                   feature.center[[4]],upstream=dd, downstream=dd,
  #                                   zeroAt=.5,type="l", 
  #                                   ylab="Averaged coverage",add=TRUE)
  # 
  # out <- cbind(out.1$density[,1],out.2$density[,2],out.3$density[,3],out.4$density[,4])
  
  #print(length(feature.center))
  #print(cvglists.ave)
  
  XX <- lapply(1:length(feature.center),function(u,feature.center,cvglists.ave,dd){
    
    #u <- 1
    sig1 <- featureAlignedSignal(cvglists.ave, feature.center[[u]], upstream=dd, downstream=dd) 
    sigs.log2 <- lapply(sig1, function(.ele) log2(.ele+1))
    out.4 <- featureAlignedDistribution2(sigs.log2, 
                                        feature.center[[u]],upstream=dd, downstream=dd,
                                        zeroAt=.5,type="l", 
                                        ylab="Averaged coverage",add=TRUE)
    
    x <- list(density=out.4$density[,u],grWidAt=out.4$grWidAt,grWidLab=out.4$grWidLab)
    x
    #out.4
  },feature.center,cvglists.ave,dd)
  
  den.l <-lapply(1:length(XX),function(u,XX){
    Z <- XX[[u]]$density
    Z
  },XX)
  out <- do.call(cbind,den.l)
  
  grWidAt.l <-lapply(1:length(XX),function(u,XX){
    Z <- XX[[u]]$grWidAt
    Z
  },XX)
  grWidAt <- do.call(cbind,grWidAt.l)
  
  grWidLab.l <-lapply(1:length(XX),function(u,XX){
    Z <- XX[[u]]$grWidLab
    Z
  },XX)
  grWidLab <- do.call(cbind,grWidLab.l)
  
  
  head(out)
  colnames(out) <- tssType
  
  drowDensity <- function(density,grWidAt,grWidLab,...) {
    dots <- list(...)
    matplot(density, ..., xaxt="n")
    axis(1, at = grWidAt, labels = grWidLab)
    lty <- if(!is.null(dots$lty)) dots$lty else 1:5
    lwd <- if(!is.null(dots$lwd)) dots$lwd else 1
    col <- if(!is.null(dots$col)) dots$col else 1:6
    legend("topright", legend=colnames(density), col=col,
           lty=lty, lwd=lwd)
  }
  
  jpeg(file.path(output.file.dir,paste0(x_name,"_density.jpeg")))
  drowDensity(out,grWidAt[,1],grWidLab[,1],type="l", lty=1,
              ylab="log2(averaged coverage)",xlab=paste0("Distance to ",class_pattern,"(bp)"),main=tfType)
  dev.off()
  
  x.pos <- seq(-dd,dd,length.out = dim(out)[1])
  out.m <- data.frame(x.pos=x.pos,out)
  
  outLong <- melt(data          = out.m,
                  id.vars       = c("x.pos"),
                  measure.vars  = tssType,
                  variable.name = "TSSByClass",
                  value.name    = "value")
  
  g <- ggplot(outLong, aes(x.pos, value,colour=TSSByClass)) + 
    geom_line() + labs(x="Distance to TSS(bp)",y="log2(averaged coverage)",title=tfType) +
    theme(plot.title=element_text(hjust=0.5)) + scale_x_continuous(name=paste0("Distance to ",class_pattern,"(bp)"), limits=c(-dd, dd))
  
  jpeg(file.path(output.file.dir,paste0(x_name,"_density_2.jpeg")))
  g
  dev.off()
  
  write.table(out.m,file = file.path(output.file.dir,paste0(x_name,"_Density_matrix_short",".txt")),
              append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T)
  
  write.table(outLong,file = file.path(output.file.dir,paste0(x_name,"_Density_matrix_long",".txt")),
              append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T)
  
  # null <- lapply(1:length(cvglists.l), function(u,cvglists.l,feature.center,dd,n,tfType){
  #   #u <- 1
  #   sig1.4.heatmap <- featureAlignedSignal(cvglists.l[[u]], feature.center[[u]], upstream=dd, downstream=dd) 
  #   sig1.4.heatmap.log2 <- lapply(sig1.4.heatmap, function(.ele) log2(.ele+1))
  #   names(sig1.4.heatmap.log2) <- paste0(tfType,"_",seq(1,n))
  #   
  #   #u <- 1
  #   #x_name_0 <- tools::file_path_sans_ext(basename(names(cvglists.l[[u]])))  
  #   #x_name_start <- str_locate_all(x_name_0,"_")[[1]][2,2] + 1
  #   #x_name_end <- str_locate_all(x_name_0,"_")[[1]][3,2] - 1
  #   
  #   x_name <- paste0(tfType,"_",names(feature.center)[u])
  #   
  #   jpeg(file.path(output.file.dir,paste0(x_name,"_",dd,"_heatmap.jpeg")))
  #   featureAlignedHeatmap(sig1.4.heatmap.log2, feature.center[[u]],sortBy=NULL,upstream=dd, downstream=dd,zeroAt=.5,res=300)
  #   dev.off()
  #   
  # },cvglists.l,feature.center,dd,n,tfType)
  
}

makeHeatMap5 <- function(){

  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  file.name.4 <- list.files(input.bed.dir,all.files = TRUE,full.names = TRUE,recursive = FALSE,include.dirs = TRUE)
  
  file.name.5 <- file.name.4[c(3,6,5,14,13)]
  
  xx <- lapply(file.name.5, function(u){
    x <- readLines(u)
    x
  })
  
  yy <- lapply(xx, function(u){
    y <- as.data.frame(do.call(rbind,strsplit(u,split="\t")),stringsAsFactors=FALSE)
    y
  })
  
  yy[[1]] <- yy[[1]][-1,]
  yy[[2]] <- yy[[2]][-140,]
  yy[[4]] <- yy[[4]][-83,]
  
  
  zz <- lapply(yy, function(u){
    if(dim(u)[2]<=3){
      z <- data.frame(u,paste0("no_",u[,1],"_",u[,2],"_",u[,3]))
    }
    else
    {
      u[grep("yes",u[,4]),4] <- paste0("yes_",u[grep("yes",u[,4]),1],"_",u[grep("yes",u[,4]),2],"_",u[grep("yes",u[,4]),3])
      
      u[-grep("yes",u[,4]),4] <- paste0("no_",u[-grep("yes",u[,4]),1],"_",u[-grep("yes",u[,4]),2],"_",u[-grep("yes",u[,4]),3])
      
      z <- u
    }
    
    colnames(z)=c("seqnames","start","end","label")
    
    z <- GRanges(z)
    z
  })
  
  data <- zz
  names(data) <- basename(file.name.5)
  
  names(data) <- gsub(".bed", "", names(data))
  
  names(data) <- gsub(".tab", "", names(data))
  
  class_patterns <- c("1_VeNR2F2repsFinal","Art_K27acUpNR2F2KD","Art_K27acDwnNR2F2KD","Vein_K27acUpNR2F2KD","Vein_K27acDwnNR2F2KD")
  file_name <- "HDAC1"
  cols= colorpanel(100,"white","lightblue","blue")
  
  #~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/OutPutCutAndRun/
    
  class_pattern <- class_patterns[1]
  sig1.4.heatmap.log2.ck.tr.Up <- generateMatrix4heatmap1(class_pattern,data,dd,input.bw.path,sortBy="Trt")
  
  
  signal.ave <- apply(sig1.4.heatmap.log2.ck.tr.Up[[1]],1,mean)
  order.index <- order(signal.ave,decreasing=T)
  sig1.4.heatmap.log2.ck.tr.Up.sorted <- sig1.4.heatmap.log2.ck.tr.Up[[1]][order.index,]
  
  
  
  ht_list <- lapply(1:length(sig1.4.heatmap.log2.ck.tr.Up), function(u,sig1.4.heatmap.log2.ck.tr.Up){
    
    sig1.4.heatmap.log2.ck.tr.1 <- sig1.4.heatmap.log2.ck.tr.Up[[u]]
  
    signal.ave <- apply(sig1.4.heatmap.log2.ck.tr.1,1,mean)
    order.index <- order(signal.ave,decreasing=T)
    sig1.4.heatmap.log2.ck.tr.1.sorted <- sig1.4.heatmap.log2.ck.tr.1[order.index,]
      
    z.mat.0 <- sig1.4.heatmap.log2.ck.tr.1.sorted
    
    #z.mat.0 <- t(scale(t(sig1.4.heatmap.log2.ck.tr.1), center=TRUE, scale=TRUE))
    
    #n.no.change <- dim(sig1.4.heatmap.log2.ck.tr.no.change[[u]])[1]
    #n.Up <- dim(sig1.4.heatmap.log2.ck.tr.Up[[u]])[1]
    #n.Down <- dim(sig1.4.heatmap.log2.ck.tr.Down[[u]])[1]
    
    ht <- Heatmap(z.mat.0, name = names( sig1.4.heatmap.log2.ck.tr.Up)[u],
                  col = cols,            
                  show_row_name = FALSE,
                  cluster_columns = F,
                  cluster_rows = F,column_names_gp = gpar(fontsize = 6))
    
    ht
    
  },sig1.4.heatmap.log2.ck.tr.Up)
  
  htList = ht_list[[1]] + ht_list[[2]] + ht_list[[3]]
  
  gb = grid.grabExpr(draw(htList))
  g <- plot_grid(gb)
  save_plot(file.path(output.file.dir,paste0(class_pattern,"_heatmap_sorted.png")),g,base_width = 10,base_height = 40)
  
  #heatmap4UpNoChangeDown(class_patterns,file_name,data, dd, input.bw.path, output.file.dir,cols=cols)
  
# 2. Heat maps
# - 3 separate heat maps at 3000 nucleotides flanking center of feature. For heat maps with more than one file, place them on the heat map in the order shown.
# 
# map HDAC1 cut and run onto (single heat map for each replicate sample):
#   
#   1_VeNR2F2repsFinal.bed
# (color of max = blue)
# 
# Art_K27acUpNR2F2KD.tab
# Art_K27acDwnNR2F2KD.tab
# (color of max=red)
# 
# Vein_K27acUpNR2F2KD.tab
# Vein_K27acDwnNR2F2KD.tab
# (color of max=blue)
  
}

# dd <- 10000
# tfType <- "HDAC1"
# input.bed.dir <- "~/Aimin/project/umw_nathan_lawson/deep_seq_data/CutAndRun/BedFilesForPlotsAndHeatMaps"
# input.bw.path <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/OutPutCutAndRun"
# class_pattern <- "TSS"
# tssType <- c("Artery","Common","Vein")
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/Results_4_CutAndRun"
# null <- overLapWithOtherFeatures5(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

# class_pattern <- "All"
# null <- overLapWithOtherFeatures5(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

# class_pattern <- "Enh"
# null <- overLapWithOtherFeatures5(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

featureAlignedDistribution2 <- function(cvglists, feature.gr, 
                                        upstream, downstream, 
                                        n.tile=100, zeroAt, ...){
  stopifnot(is(feature.gr, "GRanges"))
  dots <- list(...)
  
  grWidr <- unique(width(feature.gr))
  if(missing(upstream) || missing(downstream)){
    if(length(grWidr)!=1){
      stop("The width of feature.gr is not identical.")
    }
    if(missing(zeroAt)) {
      zeroAt <- 0.5
      message("zero is set as the center of the feature.gr")
    }
  }else{
    if(!is.numeric(upstream) || !is.numeric(downstream)){
      stop("upstream and downstream must be integers")
    }
    if(upstream<0 || downstream<0){
      stop("upstream and downstream must be not less than 0")
    }
    upstream <- as.integer(upstream)
    downstream <- as.integer(downstream)
    if(length(grWidr)!=1 || any(grWidr!=1)){
      start(feature.gr) <- start(feature.gr)+floor(width(feature.gr)/2)
      width(feature.gr) <- 1
      warning("feature.gr is set to the center of feature.gr")
    }
    if(!missing(zeroAt)){
      warning("zeroAt will be ignored.")
    }
    zeroAt <- upstream/(upstream + downstream)
    end(feature.gr) <- start(feature.gr) + downstream
    start(feature.gr) <- start(feature.gr) - upstream
    grWidr <- unique(width(feature.gr))
  }
  stopifnot(is.numeric(zeroAt))
  stopifnot(zeroAt>=0)
  if(zeroAt<=1){
    zero <- round(grWidr*zeroAt)
  }else{
    zero <- round(zeroAt)
  }
  
  grWid <- c(0, grWidr) - zero
  grWidLab <- grid.pretty(grWid)
  grWidAt <- (grWidLab+zero)/grWidr*n.tile
  if(inherits(cvglists, c("SimpleRleList", "RleList", "CompressedRleList"))){
    cvglistsName <- substitute(deparse(cvglists))
    cvglists <- list(cvglists)
    names(cvglists) <- cvglistsName
  }
  if(!is(cvglists, "list")){
    stop("cvglists must be output of featureAlignedSignal or 
         a list of SimpleRleList or RleList")
  }
  cls <- sapply(cvglists, class)
  if(all(cls=="matrix")){
    cov <- cvglists
    if(ncol(cov[[1]])!=n.tile){
      stop("n.tile must keep same as featureAligendSignal.")
    }
  }else{
    if(any(!cls %in% c("SimpleRleList", "RleList", "CompressedRleList")))
      stop("cvglists must be a list of SimpleRleList or RleList")
    cov <- featureAlignedSignal(cvglists, feature.gr, n.tile=n.tile)
  }
  ## normalized read density
  if(any(sapply(cov, function(.ele) any(is.na(.ele))))){
    warning("cvglists contain NA values. ", 
            "NA value will be omit.")
  }
  density <- sapply(cov, colMeans, na.rm=TRUE)
  try({
    matplot(density, ..., xaxt="n")
    axis(1, at = grWidAt, labels = grWidLab)
    lty <- if(!is.null(dots$lty)) dots$lty else 1:5
    lwd <- if(!is.null(dots$lwd)) dots$lwd else 1
    col <- if(!is.null(dots$col)) dots$col else 1:6
    legend("topright", legend=colnames(density), col=col,
           lty=lty, lwd=lwd)
  })
  
  res <- list(density=density,grWidAt=grWidAt,grWidLab=grWidLab)
  res
  #return(invisible(density))
  }

generateMatrix4heatmap1 <- function(class_pattern,data,dd,input.bw.path,sortBy) {
  
  tss <- names(data)[grep(class_pattern,names(data),ignore.case = T)]
  print(tss)
  
  tss.bed <- data[which(names(data) %in% tss)]
  dd <- 3000
  feature.center <- lapply(1:length(tss.bed), function(u,tss.bed){
    x <- tss.bed[[u]]
    x
  },tss.bed)
  
  names(feature.center) <- names(tss.bed)
  
  # get the regions of around feature middle +/- dd 
  feature.recentered.l <- lapply(1:length(tss.bed), function(u,tss.bed,dd){
    
    features <- tss.bed[[u]]
    
    wid <- width(features)
    feature.recentered <- feature.center <- features
    
    start(feature.center) <- start(features) + floor(wid/2)
    
    width(feature.center) <- 1
    
    start(feature.recentered) <- start(feature.center) - dd
    end(feature.recentered) <- end(feature.center) + dd
    
    feature.recentered
    
  },tss.bed,dd)
  names(feature.recentered.l) <- names(tss.bed)
  
  cvglists.l <- lapply(1:length(feature.recentered.l), function(u,feature.recentered.l,input.bw.path){
    
    files <- dir(input.bw.path, "bw")
    if(.Platform$OS.type != "windows"){
      cvglists <- sapply(file.path(input.bw.path,files), import, 
                         format="BigWig", 
                         which=feature.recentered.l[[u]], 
                         as="RleList")
    }else{## rtracklayer can not import bigWig files on Windows
      load(file.path(input.bw.path, "cvglist.rds"))
    }
    
  },feature.recentered.l,input.bw.path)
  
  tfType <- "HDAC1"
  cvglists.ck <- lapply(1:length(cvglists.l), function(u,cvglists.l){
    z <- cvglists.l[[u]][grep(tfType,names(cvglists.l[[u]]))]
    z
  },cvglists.l)
  cvglists.ul.ck <- unlist(cvglists.ck)
  
  print(names(cvglists.ul.ck))
  
  n <- length(cvglists.ul.ck)
  
  sig1.4.heatmap.ck <- featureAlignedSignal2(cvglists.ck[[1]], feature.center[[1]], upstream=dd, downstream=dd)
  sig1.4.heatmap.log2.ck <- lapply(sig1.4.heatmap.ck, function(.ele) log2(.ele+1))
  names(sig1.4.heatmap.log2.ck) <- paste0(tfType,"_",seq(1,n))
  print(names(sig1.4.heatmap.log2.ck))
  
  # tfType <- "shNR2F2"
  # cvglists.tr <- lapply(1:length(cvglists.l), function(u,cvglists.l){
  #   z <- cvglists.l[[u]][grep(tfType,names(cvglists.l[[u]]))]
  #   z
  # },cvglists.l)
  # cvglists.ul.tr <- unlist(cvglists.tr)
  # print(names(cvglists.ul.tr))
  # n <- length(cvglists.ul.tr)
  # 
  # sig1.4.heatmap.tr <- featureAlignedSignal2(cvglists.tr[[1]], feature.center[[1]], upstream=dd, downstream=dd)
  # sig1.4.heatmap.log2.tr <- lapply(sig1.4.heatmap.tr, function(.ele) log2(.ele+1))
  # names(sig1.4.heatmap.log2.tr) <- paste0(tfType,"_",seq(1,n))
  # print(names(sig1.4.heatmap.log2.tr))
  # 
  # dim(sig1.4.heatmap.log2.tr[[1]])
  # 
  # #sortBy = "Trt"
  # sig1.4.heatmap.log2.ck.tr.l <- lapply(1:length(sig1.4.heatmap.log2.tr),function(u,sig1.4.heatmap.log2.ck,sig1.4.heatmap.log2.tr,sortBy){
  #   if(sortBy=="Trt"){
  #     density.sum <- apply(sig1.4.heatmap.log2.tr[[u]],1,sum)
  #     sig1.4.heatmap.log2.tr.sorted<- sig1.4.heatmap.log2.tr[[u]][order(density.sum,decreasing = T),]        
  #     sig1.4.heatmap.log2.ck.sorted <- sig1.4.heatmap.log2.ck[[u]][match(rownames(sig1.4.heatmap.log2.tr.sorted),rownames(sig1.4.heatmap.log2.ck[[u]])),]
  #     sig1.4.heatmap.log2.ck.tr <- cbind(sig1.4.heatmap.log2.ck.sorted,sig1.4.heatmap.log2.tr.sorted)
  #   }else{
  #     density.sum <- apply(sig1.4.heatmap.log2.ck[[u]],1,sum)
  #     sig1.4.heatmap.log2.ck.sorted<- sig1.4.heatmap.log2.ck[[u]][order(density.sum,decreasing = T),]        
  #     sig1.4.heatmap.log2.tr.sorted <- sig1.4.heatmap.log2.tr[[u]][match(rownames(sig1.4.heatmap.log2.ck.sorted),rownames(sig1.4.heatmap.log2.tr[[u]])),]
  #     sig1.4.heatmap.log2.ck.tr <- cbind(sig1.4.heatmap.log2.ck.sorted,sig1.4.heatmap.log2.tr.sorted)
  #   }
  #   sig1.4.heatmap.log2.ck.tr
  # },sig1.4.heatmap.log2.ck,sig1.4.heatmap.log2.tr,sortBy)
  # 
  # sig1.4.heatmap.log2.ck.tr.l
  sig1.4.heatmap.log2.ck
}

heatmap4UpDown <- function(class_patterns,file_name,data, dd, input.bw.path, output.file.dir,selected.genes=NULL,cols=colorpanel(30,"white","white","red")) {
  
  class_pattern <- class_patterns[1]
  sig1.4.heatmap.log2.ck.tr.Up <- generateMatrix4heatmap(class_pattern,data,dd,input.bw.path,sortBy="Trt")
  
  lbsUp <- as.character(mcols(data[[grep(class_pattern,names(data))]])$label)
  
  #print(names(sig1.4.heatmap.log2.ck.tr.Up))
  
  class_pattern <- class_patterns[2]
  sig1.4.heatmap.log2.ck.tr.Down <- generateMatrix4heatmap(class_pattern,data,dd,input.bw.path,sortBy="CK")
  
  lbsDown <- as.character(mcols(data[[grep(class_pattern,names(data))]])$label)
  #data[[grep(class_pattern,names(data))]]
  
  lbs <- c(lbsUp,lbsDown)
  
  if(is.null(selected.genes)) {selected.genes <- lbs[grep("yes",lbs)]}
  
  ht_list <- lapply(1:length(sig1.4.heatmap.log2.ck.tr.Up), function(u,sig1.4.heatmap.log2.ck.tr.Down,selected.genes){
    
    sig1.4.heatmap.log2.ck.tr.1 <- rbind(sig1.4.heatmap.log2.ck.tr.Up[[u]],sig1.4.heatmap.log2.ck.tr.Down[[u]])
  
    z.mat.0 <- t(scale(t(sig1.4.heatmap.log2.ck.tr.1), center=TRUE, scale=TRUE))
    
    print(rownames(z.mat.0))
    
    labels <- as.character(rownames(z.mat.0))
    selected.index <- match(selected.genes,labels)
    markers <- labels[selected.index]
    
    marker_idx <- selected.index
    
    print(labels[marker_idx])
    
    n.Up <- dim(sig1.4.heatmap.log2.ck.tr.Up[[u]])[1]
    n.Down <- dim(sig1.4.heatmap.log2.ck.tr.Down[[u]])[1]
    
    ht <- Heatmap(z.mat.0, name = paste0("z-score","-rep",u),
                  #col = rev(redblue(30))[-seq(35, 35)],  
                  #col = rev(redblue(30)),
                  # col = colorpanel(30,"white","white","red"),
                  col = cols,
                  show_row_name = FALSE,
                  cluster_columns = F,
                  cluster_rows = F,column_names_gp = gpar(fontsize = 6),column_split = rep(c("Control", "NR2F2KD"),c(100,100)),row_split = factor(rep(c("Up","Down"),c(n.Up,n.Down)),levels = c("Up","Down"))) + rowAnnotation(link = anno_mark(at = marker_idx, labels = markers),width = unit(5, "mm") +        max_text_width(markers))
    
    ht
    
  },sig1.4.heatmap.log2.ck.tr.Down,selected.genes)
  
  htList = ht_list[[1]] + ht_list[[2]] + ht_list[[3]]
  
  gb = grid.grabExpr(draw(htList))
  g <- plot_grid(gb)
  save_plot(file.path(output.file.dir,paste0(file_name,"_heatmap.png")),g,base_width = 25,base_height = 40)
}

generateData4DensityPlot <- function(class_pattern,data,dd,input.bw.path,sortBy) {
  
  tss <- names(data)[grep(class_pattern,names(data),ignore.case = T)]
  print(tss)
  
  tss.bed <- data[which(names(data) %in% tss)]
  
  dd <- 3000
  feature.center <- lapply(1:length(tss.bed), function(u,tss.bed){
    x <- tss.bed[[u]]
    x
  },tss.bed)
  
  names(feature.center) <- names(tss.bed)
  
  # get the regions of around feature middle +/- dd 
  feature.recentered.l <- lapply(1:length(tss.bed), function(u,tss.bed,dd){
    
    features <- tss.bed[[u]]
    
    wid <- width(features)
    feature.recentered <- feature.center <- features
    
    start(feature.center) <- start(features) + floor(wid/2)
    
    width(feature.center) <- 1
    
    start(feature.recentered) <- start(feature.center) - dd
    end(feature.recentered) <- end(feature.center) + dd
    
    feature.recentered
    
  },tss.bed,dd)
  names(feature.recentered.l) <- names(tss.bed)
  
  cvglists.l <- lapply(1:length(feature.recentered.l), function(u,feature.recentered.l,input.bw.path){
    
    files <- dir(input.bw.path, "bw")
    if(.Platform$OS.type != "windows"){
      cvglists <- sapply(file.path(input.bw.path,files), import, 
                         format="BigWig", 
                         which=feature.recentered.l[[u]], 
                         as="RleList")
    }else{## rtracklayer can not import bigWig files on Windows
      load(file.path(input.bw.path, "cvglist.rds"))
    }
    
  },feature.recentered.l,input.bw.path)
  
  cvglists.ul <- unlist(cvglists.l)
  
  n <- length(cvglists.ul)/length(tss.bed)
  
  if(n==2){
    cvglists.ave <- lapply(1:length(tss.bed),function(u,cvglists.ul){
      
      Class1 <- (cvglists.ul[[2*u-1]]+cvglists.ul[[2*u]])/n
      
      Class1
      
    },cvglists.ul)
  }
  
  if(n==3){
    cvglists.ave <- lapply(1:length(tss.bed),function(u,cvglists.ul){
      
      Class1 <- (cvglists.ul[[3*u-2]]+cvglists.ul[[3*u-1]]+cvglists.ul[[3*u]])/n
      
      Class1
      
    },cvglists.ul)
  }
  
  names(cvglists.ave) <- names(tss.bed)
  
  XX <- lapply(1:length(feature.center),function(u,feature.center,cvglists.ave,dd){
    
    sig1 <- featureAlignedSignal(cvglists.ave, feature.center[[u]], upstream=dd, downstream=dd) 
    sigs.log2 <- lapply(sig1, function(.ele) log2(.ele+1))
    out.4 <- featureAlignedDistribution2(sigs.log2, 
                                         feature.center[[u]],upstream=dd, downstream=dd,
                                         zeroAt=.5,type="l", 
                                         ylab="Averaged coverage",add=TRUE)
    
    x <- list(density=out.4$density[,u],grWidAt=out.4$grWidAt,grWidLab=out.4$grWidLab)
    x

  },feature.center,cvglists.ave,dd)
  
  den.l <-lapply(1:length(XX),function(u,XX){
    Z <- XX[[u]]$density
    Z
  },XX)
  out <- do.call(cbind,den.l)
  
  grWidAt.l <-lapply(1:length(XX),function(u,XX){
    Z <- XX[[u]]$grWidAt
    Z
  },XX)
  grWidAt <- do.call(cbind,grWidAt.l)
  
  grWidLab.l <-lapply(1:length(XX),function(u,XX){
    Z <- XX[[u]]$grWidLab
    Z
  },XX)
  grWidLab <- do.call(cbind,grWidLab.l)
  
  #colnames(out) <- names(tss.bed)
  
  colnames(out) <- "VeNR2F2repsFinal"
  
  drawDensity <- function(density,grWidAt,grWidLab,...) {
    dots <- list(...)
    matplot(density, ..., xaxt="n")
    axis(1, at = grWidAt, labels = grWidLab)
    lty <- if(!is.null(dots$lty)) dots$lty else 1:5
    lwd <- if(!is.null(dots$lwd)) dots$lwd else 1
    col <- if(!is.null(dots$col)) dots$col else 1:6
    legend("topright", legend=colnames(density), col=col,
           lty=lty, lwd=lwd)
  }
  
  x_name <- paste0(tfType,"_",dd,"_",names(tss.bed))
  
  jpeg(file.path(output.file.dir,paste0(x_name,"_density.jpeg")))
  drawDensity(out,grWidAt[,1],grWidLab[,1],type="l", lty=1,
              ylab="log2(averaged coverage)",xlab=paste0("Distance to ",class_pattern,"(bp)"),main=tfType)
  dev.off()
  
  x.pos <- seq(-dd,dd,length.out = dim(out)[1])
  out.m <- data.frame(x.pos=x.pos,out)
  tssType <- "VeNR2F2repsFinal"
  outLong <- melt(data          = out.m,
                  id.vars       = c("x.pos"),
                  measure.vars  = tssType,
                  variable.name = "TSSByClass",
                  value.name    = "value")
  
  g <- ggplot(outLong, aes(x.pos, value,colour=TSSByClass)) + 
    geom_line() + labs(x="Distance to TSS(bp)",y="log2(averaged coverage)",title=tfType) +
    theme(plot.title=element_text(hjust=0.5)) + scale_x_continuous(name=paste0("Distance to ",class_pattern,"(bp)"), limits=c(-dd, dd))
  
  jpeg(file.path(output.file.dir,paste0(x_name,"_density_2.jpeg")))
  g
  dev.off()
  
  write.table(out.m,file = file.path(output.file.dir,paste0(x_name,"_Density_matrix_short",".txt")),append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T)
  
  write.table(outLong,file = file.path(output.file.dir,paste0(x_name,"_Density_matrix_long",".txt")),append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T)
  
}

generateMatrix4heatmap2 <- function(class_pattern,data,dd,input.bw.path,sortBy) {
  
  tss <- names(data)[grep(class_pattern,names(data),ignore.case = T)]
  print(tss)
  
  tss.bed <- data[which(names(data) %in% tss)]
  
  feature.center <- lapply(1:length(tss.bed), function(u,tss.bed,dd){
    x <- tss.bed[[u]]
    x
  },tss.bed,dd)
  
  names(feature.center) <- names(tss.bed)
  
  # get the regions of around feature middle +/- dd 
  feature.recentered.l <- lapply(1:length(tss.bed), function(u,tss.bed,dd){
    
    features <- tss.bed[[u]]
    
    wid <- width(features)
    feature.recentered <- feature.center <- features
    
    start(feature.center) <- start(features) + floor(wid/2)
    
    width(feature.center) <- 1
    
    start(feature.recentered) <- start(feature.center) - dd
    end(feature.recentered) <- end(feature.center) + dd
    
    feature.recentered
    
  },tss.bed,dd)
  names(feature.recentered.l) <- names(tss.bed)
  
  cvglists.l <- lapply(1:length(feature.recentered.l), function(u,feature.recentered.l,input.bw.path){
    
    files <- dir(input.bw.path, "bigwig")
    if(.Platform$OS.type != "windows"){
      cvglists <- sapply(file.path(input.bw.path,files), import, 
                         format="BigWig", 
                         which=feature.recentered.l[[u]], 
                         as="RleList")
    }else{## rtracklayer can not import bigWig files on Windows
      load(file.path(input.bw.path, "cvglist.rds"))
    }
    
  },feature.recentered.l,input.bw.path)
  
  tfType <- "V29NR2F2in_"
  cvglists.ck <- lapply(1:length(cvglists.l), function(u,cvglists.l){
    z <- cvglists.l[[u]][grep(tfType,names(cvglists.l[[u]]))]
    z
  },cvglists.l)
  cvglists.ul.ck <- unlist(cvglists.ck)
  
  print(names(cvglists.ul.ck))
  
  n <- length(cvglists.ul.ck)
  
  sig1.4.heatmap.ck <- featureAlignedSignal(cvglists.ck[[1]], feature.center[[1]], upstream=dd, downstream=dd)
  sig1.4.heatmap.log2.ck <- lapply(sig1.4.heatmap.ck, function(.ele) log2(.ele+1))
  names(sig1.4.heatmap.log2.ck) <- paste0(tfType,"_",seq(1,n))
  print(names(sig1.4.heatmap.log2.ck))
  
  tfType <- "V29NR2F2_"
  cvglists.tr <- lapply(1:length(cvglists.l), function(u,cvglists.l){
    z <- cvglists.l[[u]][grep(tfType,names(cvglists.l[[u]]))]
    z
  },cvglists.l)
  cvglists.ul.tr <- unlist(cvglists.tr)
  print(names(cvglists.ul.tr))
  n <- length(cvglists.ul.tr)
  
  sig1.4.heatmap.tr <- featureAlignedSignal(cvglists.tr[[1]], feature.center[[1]], upstream=dd, downstream=dd)
  sig1.4.heatmap.log2.tr <- lapply(sig1.4.heatmap.tr, function(.ele) log2(.ele+1))
  names(sig1.4.heatmap.log2.tr) <- paste0(tfType,"_",seq(1,n))
  print(names(sig1.4.heatmap.log2.tr))
  
  dim(sig1.4.heatmap.log2.tr[[1]])
  
  #sortBy = "Trt"
  sig1.4.heatmap.log2.ck.tr.l <- lapply(1:length(sig1.4.heatmap.log2.tr),function(u,sig1.4.heatmap.log2.ck,sig1.4.heatmap.log2.tr,sortBy){
    if(sortBy=="Trt"){
      density.sum <- apply(sig1.4.heatmap.log2.tr[[u]],1,sum)
      sig1.4.heatmap.log2.tr.sorted<- sig1.4.heatmap.log2.tr[[u]][order(density.sum,decreasing = T),]        
      sig1.4.heatmap.log2.ck.sorted <- sig1.4.heatmap.log2.ck[[u]][match(rownames(sig1.4.heatmap.log2.tr.sorted),rownames(sig1.4.heatmap.log2.ck[[u]])),]
      sig1.4.heatmap.log2.ck.tr <- cbind(sig1.4.heatmap.log2.ck.sorted,sig1.4.heatmap.log2.tr.sorted)
    }else{
      density.sum <- apply(sig1.4.heatmap.log2.ck[[u]],1,sum)
      sig1.4.heatmap.log2.ck.sorted<- sig1.4.heatmap.log2.ck[[u]][order(density.sum,decreasing = T),]        
      sig1.4.heatmap.log2.tr.sorted <- sig1.4.heatmap.log2.tr[[u]][match(rownames(sig1.4.heatmap.log2.ck.sorted),rownames(sig1.4.heatmap.log2.tr[[u]])),]
      sig1.4.heatmap.log2.ck.tr <- cbind(sig1.4.heatmap.log2.ck.sorted,sig1.4.heatmap.log2.tr.sorted)
    }
    sig1.4.heatmap.log2.ck.tr
  },sig1.4.heatmap.log2.ck,sig1.4.heatmap.log2.tr,sortBy)
  
  sig1.4.heatmap.log2.ck.tr.l
  
}

drawMat4HeatMap <- function(sig1.4.heatmap.log2.ck.tr.Up, cols) {
  ht_list <- lapply(1:length(sig1.4.heatmap.log2.ck.tr.Up), function(u,cols){
    
    #u <- 1
    sig1.4.heatmap.log2.ck.tr.1 <- sig1.4.heatmap.log2.ck.tr.Up[[u]]
    
  
    #z.mat.0 <- t(scale(t(sig1.4.heatmap.log2.ck.tr.1), center=TRUE, scale=TRUE))
    
    z.mat.0 <- sig1.4.heatmap.log2.ck.tr.1
    
    print(rownames(z.mat.0))
    
    ht <- Heatmap(z.mat.0, name = paste0("z-score","-rep",u),
                  #col = rev(redblue(30))[-seq(35, 35)],  
                  #col = rev(redblue(30)),
                  # col = colorpanel(30,"white","white","red"),
                  col = cols,
                  show_row_name = FALSE,
                  cluster_columns = F,
                  cluster_rows = F,column_names_gp = gpar(fontsize = 6),column_split = rep(c(paste0("V29NR2F2in_",u),paste0("V29NR2F2_",u)),c(100,100)))
    
    ht
    
  },cols)
  
  htList = ht_list[[1]] + ht_list[[2]]
  htList
}

makeGrDataFromBedFiles <- function(input.bed.dir) {

  # ArteryAll.bed, CommonAll.bed,VeinAll.bed
  # ArteryTSS.bed,CommonTSS.bed, VeinTSS.bed
  # ArteryEnh.bed,CommonEnh.bed, VeinEnh.bed
  
  file.name.4 <- list.files(input.bed.dir,pattern=".bed$",all.files = TRUE,full.names = TRUE,recursive = FALSE,include.dirs = TRUE)
  
  xx <- lapply(file.name.4, function(u){
    x <- readLines(u)
    x
  })
  
  yy <- lapply(xx, function(u){
    y <- as.data.frame(do.call(rbind,strsplit(u,split="\t")),stringsAsFactors=FALSE)
    y
  })
  
  yy[[4]] <- yy[[4]][-1,]
  yy[[7]] <- yy[[7]][-1,]
  yy[[10]] <- yy[[10]][-1,]
  
  yy[[3]] <- yy[[3]][-1,]
  yy[[6]] <- yy[[6]][-1,]
  yy[[9]] <- yy[[9]][-1,]
  
  names(yy) <- tools::file_path_sans_ext(basename(file.name.4))
  
  zz <- lapply(yy[c(2,5,8,4,7,10,3,6,9)], function(u){
    z <- data.frame(u[,c(1:3)],paste0(as.character(u[,1]), ":",u[,2], "-",u[,3]))
    colnames(z)=c("seqnames","start","end","label")
    z <- GRanges(z)
    z
  })
  zz
}

input.bed.dir <- "~/Aimin/project/umw_nathan_lawson/deep_seq_data/CutAndRun/BedFilesForPlotsAndHeatMaps"
data <- makeGrDataFromBedFiles(input.bed.dir)
class_patterns_1 <- names(data)
output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/Results_4_CutAndRun/V29NR2F2"
input.bw.path <- "~/Aimin/nl/umw_nathan_lawson/NR2F2bigwigs"

if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}

lapply(1:length(class_patterns_1), function(u,class_patterns_1){
  
  class_pattern <- class_patterns_1[u]
  
  NR2F2bigwigs.data <- generateMatrix4heatmap2(class_pattern,data,dd,input.bw.path,sortBy)
  
  V29NR2F2.ht <- drawMat4HeatMap(NR2F2bigwigs.data, cols)
  
  gb = grid.grabExpr(draw(V29NR2F2.ht))
  g <- plot_grid(gb)
  save_plot(file.path(output.file.dir,paste0(class_pattern,"_heatmap.png")),g,base_width = 25,base_height = 40)
  
},class_patterns_1)

save.image(file = file.path(output.file.dir,"cut_run_data.RData"))
savehistory(file = file.path(output.file.dir,"cut_run_data.Rhistory"))

save.image(file = file.path(output.file.dir,"NR2F2bigwigs.RData"))
savehistory(file = file.path(output.file.dir,"NR2F2bigwigs.Rhistory"))

