#!/usr/bin/env Rscript

# "Dropbox\pub_old_galaxy\F121-9\parsed NADs\129 adjPlt0.05-1.5FC-countF200-minimump.xls"

bed.129 <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/parsed NADs/129_adjPlt0.05-1.5FC-countF200-minimump.txt"

# "Dropbox\pub_old_galaxy\F121-9\parsed NADs\Cast adjPlt0.05-1.5FC-countF200-minimump.xls"

bed.cast <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/parsed NADs/Cast_adjPlt0.05-1.5FC-countF200-minimump.txt"

bed.files <- c(bed.129,bed.cast)

bed.in<-lapply(1:length(bed.files),function(u,bed.files){
  x <- bed.files[[u]]
  peaks=read.table(x,header=T)
  colnames(peaks)[1:3]= c("chr","start","end")
  peaks=toGRanges(peaks)
  peaks
},bed.files)

names(bed.in) <- c("F129","Cast")

output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/DiffBindAnalysis/peaks"

outGrl(bed.in,speci="Mm",output.file.dir)
  
bam.dir <- "~/Aimin/DropboxUmass/pub_old_galaxy/F121-9/Aizhan.allele.specific/DNASEQ/hybrid_align/markDup_align"

file.1 <- list.files(bam.dir,pattern="markDup.bam",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)

index.genome1 <- grep("genome1",basename(file.1))
index.genome2 <- grep("genome2",basename(file.1))

bam.genome1 <- file.1[index.genome1][c(1,3,17,19)]
bam.genome2 <- file.1[index.genome2][c(1,3,17,19)]

tamoxifen <- dba(sampleSheet="tamoxifen.csv", dir=basedir,peakCaller="macs", peakFormat="raw", scoreCol=5 )

# prepare sample.csv for dba in DiffBind

basedir1 <- "~/Aimin/DropboxUmass/NADfinder/Aimin/DiffBindAnalysis"

sample.csv.1 <- data.frame(SampleID=c("Cast_1","Cast_2","F129_1","F129_2"),	Condition=c(rep("Cast",2),rep("F129",2)),Replicate=rep(c("1","2"),2),bamReads=c(file.path("reads",basename(bam.genome1[3])),file.path("reads",basename(bam.genome1[4])),file.path("reads",basename(bam.genome2[3])),file.path("reads",basename(bam.genome2[4]))),ControlID=c("g1gDNA1","g1gDNA2","g2gDNA1","g2gDNA2"),bamControl=c(file.path("reads",basename(bam.genome1[1])),file.path("reads",basename(bam.genome1[2])),file.path("reads",basename(bam.genome2[1])),file.path("reads",basename(bam.genome2[2]))),	Peaks=c(rep(file.path("peaks","Cast.bed"),2),rep(file.path("peaks","F129.bed"),2)),PeakCaller=rep("NAD",4))

bam.genome1
bam.genome2

write.table(sample.csv.1,file = file.path(basedir1,paste0("sample",".csv")),
            append = FALSE, quote = F, sep = ",",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T)

tamoxifen <- dba(sampleSheet="sample.csv", dir=basedir1,peakCaller="NAD")
tamoxifen.paul.25k <- dba.count(tamoxifen, summits=25000)
tamoxifen.paul.25k.1 <- dba.contrast(tamoxifen.paul.25k, categories=DBA_CONDITION,minMembers=2)
tamoxifen.paul.25k.2 <- dba.analyze(tamoxifen.paul.25k.1)
tamoxifen.DB.25k.paul <- dba.report(tamoxifen.paul.25k.2,th=1, bCalled=TRUE)

dba.plotMA(tamoxifen.paul.25k.2)
dba.plotVolcano(tamoxifen.paul.25k.2)
pvals <- dba.plotBox(tamoxifen.paul.25k.2)

corvals <- dba.plotHeatmap(tamoxifen.paul.25k.2, contrast=1, correlations=FALSE)

write.table(tamoxifen.DB.25k.paul,file = file.path(basedir1,paste0("DiffBind_25K",".csv")),
            append = FALSE, quote = F, sep = ",",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T)

tamoxifen.paul <- dba.count(tamoxifen, summits=250)
tamoxifen.paul.50k <- dba.count(tamoxifen, summits=50000)
tamoxifen.paul.50k.1 <- dba.contrast(tamoxifen.paul.50k, categories=DBA_CONDITION,minMembers=2)
tamoxifen.paul.50k.2 <- dba.analyze(tamoxifen.paul.50k.1)
tamoxifen.DB.50k.paul <- dba.report(tamoxifen.paul.50k.2,th=1, bCalled=TRUE)

dba.plotMA(tamoxifen.paul.50k.2)
write.table(tamoxifen.DB.50k.paul,file = file.path(basedir1,paste0("DiffBind_50K",".csv")),
            append = FALSE, quote = F, sep = ",",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T)

plot(tamoxifen)
tamoxifen.paul.1 <- dba.contrast(tamoxifen.paul, categories=DBA_CONDITION,minMembers=2)

tamoxifen.paul.2 <- dba.analyze(tamoxifen.paul.1)

tamoxifen.DB.paul <- dba.report(tamoxifen.paul.2,th=1, bCalled=TRUE)

write.table(tamoxifen.DB.paul,file = file.path(basedir1,paste0("DiffBind",".csv")),
            append = FALSE, quote = F, sep = ",",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T)

dba.plotPCA(tamoxifen.paul.1,label=DBA_CONDITION)
dba.plotPCA(tamoxifen.paul.2, contrast=1,label=DBA_CONDITION)
dba.plotMA(tamoxifen.paul.2)

example <- data(tamoxifen_counts)
tamoxifen <- dba.contrast(tamoxifen, categories=DBA_CONDITION)

tamoxifen <- dba.analyze(tamoxifen)
tamoxifen
tamoxifen.DB <- dba.report(tamoxifen)

save.image(file = file.path(dirname(output.file.dir),paste0(basename(output.file.dir),".RData")))
savehistory(file = file.path(dirname(output.file.dir),paste0(basename(output.file.dir),".Rhistory")))
