#!/usr/bin/env Rscript

# For the data of Shan
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("HTqPCR", version = "3.8")

library(HTqPCR)
library(readxl)

all.R.commands <- system.file("doc", "HTqPCR.Rnw",package = "HTqPCR")
Stangle(all.R.commands)

input.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/Shan/"
output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/Shan/Output"

if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}

fastq.files <- list.files(path = input.file.dir,pattern="*xlsx$",full.names = TRUE, recursive = TRUE)

XX <- lapply(fastq.files[c(1,2)], function(u){
  x <- read_excel(u)
  x
})

XXX <- XX[[1]]

tr <- unique(XXX$Treatment)
cd <- unique(XXX$CD11c)
sp <- unique(XXX$species)
an <- unique(XXX$animal)
fl <- unique(XXX$FLUIDIGM)

trSpCombination <- function(XXX,cd,trP,spP) {
  
  temp <- XXX[which(XXX$Treatment==trP&XXX$species==spP),]
  n <- dim(temp)[2]
  
  if(length(unique(temp$CD11c))==2)
  {
  re <- lapply(6:n, function(u,temp,cd){
    
  x <- as.numeric(as.data.frame(temp[temp$CD11c==cd[1],])[,u])
  y <- as.numeric(as.data.frame(temp[temp$CD11c==cd[2],])[,u])
  
  if(var(x)==0&var(y)==0){
    cv1  <- var(x)/mean(x)
    cv2 <- var(y)/mean(y)
    rd <- data.frame(gene=colnames(temp)[u],x=mean(x),y=mean(y),pValue=NA)
  }else{
    stat <- t.test(x,y)
    rd <- data.frame(gene=colnames(temp)[u],x=as.numeric(stat$estimate[1]),y=as.numeric(stat$estimate[2]),pValue=stat$p.value)
  }
  rd
  
  },temp,cd)
  
  #re2 <- re[-which(sapply(re, is.null))]
  re3 <- do.call(rbind,re)
  colnames(re3) <- c("Gene",cd[1],cd[2],"pValue")
  trP <- gsub("/","-",trP)
  filename <- paste0(cd[1],"-vs-",cd[2],"-at-",trP,"-",spP,".csv")
  
  pdf(file=file.path(output.file.dir,paste0(cd[1],"-vs-",cd[2],"-at-",trP,"-",spP,".pdf")))
  hist(re3$pValue)
  dev.off()
  
  write.table(re3,file=file.path(output.file.dir,filename),sep=",",quote = FALSE,row.names = FALSE,col.names = TRUE)
  }
  
}

trSpCombination(XXX,cd,tr[1],sp[1])
trSpCombination(XXX,cd,tr[1],sp[2])

trSpCombination(XXX,cd,tr[2],sp[1])
trSpCombination(XXX,cd,tr[2],sp[2])

trSpCombination(XXX,cd,tr[3],sp[1])
trSpCombination(XXX,cd,tr[3],sp[2])

trSpCombination(XXX,cd,tr[4],sp[1])
trSpCombination(XXX,cd,tr[4],sp[2])

trSpCombination(XXX,cd,tr[5],sp[1])
trSpCombination(XXX,cd,tr[5],sp[2])

null <- lapply(1:length(tr), function(u,XXX,tr,sp){
  
  s <- tr[u]
  null <- lapply(1:length(sp),function(u,XXX,s,sp){
    trSpCombination(XXX,cd,s,sp[u])
  },XXX,s,sp)

},XXX,tr,sp)

trcdCombination <- function(XXX,sp,trP,cdP) {
  
  temp <- XXX[which(XXX$Treatment==trP&XXX$CD11c==cdP),]
  n <- dim(temp)[2]
  
  if(length(unique(temp$species))==2)
  {
    re <- lapply(6:n, function(u,temp,sp){
      x <- as.numeric(as.data.frame(temp[temp$species==sp[1],])[,u])
      y <- as.numeric(as.data.frame(temp[temp$species==sp[2],])[,u])
      if(var(x)==0&var(y)==0){
        cv1  <- var(x)/mean(x)
        cv2 <- var(y)/mean(y)
        rd <- data.frame(gene=colnames(temp)[u],x=mean(x),y=mean(y),pValue=NA)
      }else{
        stat <- t.test(x,y)
        rd <- data.frame(gene=colnames(temp)[u],x=as.numeric(stat$estimate[1]),y=as.numeric(stat$estimate[2]),pValue=stat$p.value)
      }
      rd
      
    },temp,sp)
    
    #re2 <- re[-which(sapply(re, is.null))]
    re3 <- do.call(rbind,re)
    colnames(re3) <- c("Gene",sp[1],sp[2],"pValue")
    trP <- gsub("/","-",trP)
    filename <- paste0(sp[1],"-vs-",sp[2],"-at-",trP,"-",cdP,".csv")
    
    pdf(file=file.path(output.file.dir,paste0(sp[1],"-vs-",sp[2],"-at-",trP,"-",cdP,".pdf")))
    hist(re3$pValue)
    dev.off()
    
    write.table(re3,file=file.path(output.file.dir,filename),sep=",",quote = FALSE,row.names = FALSE,col.names = TRUE)
  }
  
}

null <- lapply(1:length(tr), function(u,XXX,tr,cd){
  
  s <- tr[u]
  null <- lapply(1:length(cd),function(u,XXX,s,cd){
    trcdCombination(XXX,sp,s,cd[u])
  },XXX,s,cd)
  
},XXX,tr,cd)


output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/Shan/Output4File2"
if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}

YYY <- XX[[2]]

unique(YYY$sample1)
tr <- unique(YYY$Treatment)
cd <- unique(YYY$mDC_DN_DC)
sp <- unique(YYY$Species)

trcdCombination1 <- function(YYY,sp,trP,cdP) {
  
  temp <- YYY[which(YYY$Treatment==trP&YYY$mDC_DN_DC==cdP),]
  n <- dim(temp)[2]
  
  if(length(unique(temp$Species))==2)
  {
    re <- lapply(7:n, function(u,temp,sp){
      x <- as.numeric(as.data.frame(temp[temp$Species==sp[1],])[,u])
      y <- as.numeric(as.data.frame(temp[temp$Species==sp[2],])[,u])
      if(var(x)==0&var(y)==0){
        cv1  <- var(x)/mean(x)
        cv2 <- var(y)/mean(y)
        rd <- data.frame(gene=colnames(temp)[u],x=mean(x),y=mean(y),pValue=NA)
      }else{
        stat <- t.test(x,y)
        rd <- data.frame(gene=colnames(temp)[u],x=as.numeric(stat$estimate[1]),y=as.numeric(stat$estimate[2]),pValue=stat$p.value)
      }
      rd
      
    },temp,sp)
    
    re3 <- do.call(rbind,re)
    colnames(re3) <- c("Gene",sp[1],sp[2],"pValue")
    trP <- gsub("/","-",trP)
    cdP <- gsub(" ","_",cdP)
    filename <- paste0(sp[1],"-vs-",sp[2],"-at-",trP,"-",cdP,".csv")
    
    pdf(file=file.path(output.file.dir,paste0(sp[1],"-vs-",sp[2],"-at-",trP,"-",cdP,".pdf")))
    hist(re3$pValue)
    dev.off()
    
    write.table(re3,file=file.path(output.file.dir,filename),sep=",",quote = FALSE,row.names = FALSE,col.names = TRUE)
  }
  
}

trcdCombination1(YYY,sp,tr[1],cd[1])
trcdCombination1(YYY,sp,tr[2],cd[1])
    
ZZZ <- read_excel(fastq.files[2],sheet = 2)

unique(ZZZ$sample1)
tr <- unique(ZZZ$Treatment)
cd <- unique(ZZZ$mDC_DN_DC)
sp <- unique(ZZZ$Species)

output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/Shan/Output4File2Sheet2"
if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}

trcdCombination1(ZZZ,sp,tr[1],cd[1])
trcdCombination1(ZZZ,sp,tr[2],cd[1])
trcdCombination1(ZZZ,sp,tr[3],cd[1])
trcdCombination1(ZZZ,sp,tr[4],cd[1])


ZZZ1 <- ZZZ[ZZZ$Treatment=="TLR2/4",]

YYY1 <- YYY[YYY$Treatment=="TLR2/4",]

comm <- intersect(colnames(YYY1), colnames(ZZZ1))

YYY2 <- YYY1[,match(comm,colnames(YYY1))]
ZZZ2 <- ZZZ1[,match(comm,colnames(ZZZ1))]

YYY.ZZZ <- rbind(YYY2,ZZZ2)

tr <- unique(YYY.ZZZ$Treatment)
cd <- unique(YYY.ZZZ$mDC_DN_DC)
sp <- unique(YYY.ZZZ$Species)

trSpCombination1 <- function(XXX,cd,trP,spP) {
  
  temp <- XXX[which(XXX$Treatment==trP&XXX$Species==spP),]
  n <- dim(temp)[2]
  
  if(length(unique(temp$mDC_DN_DC))==2)
  {
    re <- lapply(7:n, function(u,temp,cd){
      
      x <- as.numeric(as.data.frame(temp[temp$mDC_DN_DC==cd[1],])[,u])
      y <- as.numeric(as.data.frame(temp[temp$mDC_DN_DC==cd[2],])[,u])
      
      if(var(x)==0&var(y)==0){
        cv1  <- var(x)/mean(x)
        cv2 <- var(y)/mean(y)
        rd <- data.frame(gene=colnames(temp)[u],x=mean(x),y=mean(y),pValue=NA)
      }else{
        stat <- t.test(x,y)
        rd <- data.frame(gene=colnames(temp)[u],x=as.numeric(stat$estimate[1]),y=as.numeric(stat$estimate[2]),pValue=stat$p.value)
      }
      rd
      
    },temp,cd)
    
    #re2 <- re[-which(sapply(re, is.null))]
    re3 <- do.call(rbind,re)
    
    cd[1] <- gsub(" ","_",cd[1])
    cd[2] <- gsub(" ","_",cd[2])
    
    colnames(re3) <- c("Gene",cd[1],cd[2],"pValue")
    trP <- gsub("/","-",trP)
    filename <- paste0(cd[1],"-vs-",cd[2],"-at-",trP,"-",spP,".csv")
    
    pdf(file=file.path(output.file.dir,paste0(cd[1],"-vs-",cd[2],"-at-",trP,"-",spP,".pdf")))
    hist(re3$pValue)
    dev.off()
    
    write.table(re3,file=file.path(output.file.dir,filename),sep=",",quote = FALSE,row.names = FALSE,col.names = TRUE)
  }
  
}

output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/Shan/Output4File2Sheet1And2"
if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}

trSpCombination1(YYY.ZZZ,cd,tr[1],sp[1])
trSpCombination1(YYY.ZZZ,cd,tr[1],sp[2])
  
save.image(file = file.path(output.file.dir,"Shan.RData"))
savehistory(file = file.path(output.file.dir,"Shan.Rhistory"))
