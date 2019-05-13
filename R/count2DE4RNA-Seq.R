#!/usr/bin/env Rscript

# For the data of HardyWilliam

input.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/Output4HardyWilliamNotRmDup"
output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/Output4HardyWilliamNotRmDup/Output"

if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}

fastq.files <- list.files(path = input.file.dir,pattern= "*txt$",full.names = TRUE, recursive = TRUE)

fastq.files <- fastq.files[grep("count.txt",fastq.files)]

XX <- lapply(fastq.files[c(1)], function(u){
  
  x <- read.table(u,header = T)
  x
})

#DE.data <- read.csv(fastq.files[5])

count.data.tli <- XX[[1]]

rownames(count.data.tli) <-  count.data.tli[,1]
count.data.tli.1 <- count.data.tli[,-(1:6)]
library(edgeR)
rawdata.tli.cpms <- cpm(count.data.tli.1)

keep <- rowSums(rawdata.tli.cpms > 0.5)  >=  2
rawdata.tli.after.filter.by.cpm <- count.data.tli.1[keep,]

tli.new.col.name <- gsub(".project.umw_minggang_fang.HardyWilliam.Output4HardyWilliamNotRmDup.Alignment.","",colnames(rawdata.tli.after.filter.by.cpm))
tli.new.col.name.1 <- gsub("X","",tli.new.col.name)
tli.new.col.name.2 <- gsub("_sorted.bam","",tli.new.col.name)
colnames(rawdata.tli.after.filter.by.cpm) <- tli.new.col.name.2

# Set DOE(Design Of Experiemnt)

tli.colData <- data.frame(SampleName=colnames(rawdata.tli.after.filter.by.cpm),batach=rep(c("b1","b1","b1"),3),condition=c(rep("ref",3),rep("low",3),rep("high",3)))

tli.colData$group <- factor(paste0(tli.colData$batach,"_",tli.colData$condition))

library(DESeq2)

ddsFullCountTable<-DESeqDataSetFromMatrix(
    countData=rawdata.tli.after.filter.by.cpm,
    colData=tli.colData,
    design= ~condition)

dds <-DESeq(ddsFullCountTable)
  
res.low.ref <- results(dds,c("condition","low","ref"))
res.low.ref.with.counts <- data.frame(res.low.ref,counts(dds,normalized = F),counts(dds,normalized = T))
  
res.high.ref <- results(dds,c("condition","high","ref"))
res.high.ref.with.counts <- data.frame(res.high.ref,counts(dds,normalized = F),counts(dds,normalized = T))
  
res.high.low <- results(dds,c("condition","high","low"))
res.high.low.with.counts <- data.frame(res.high.low,counts(dds,normalized = F),counts(dds,normalized = T))
  
library(biomaRt)

mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol"),values=rownames(res.low.ref.with.counts),mart= mart)

getGeneSymbolOut <- function(res.low.ref.with.counts,G_list,fileName) {
  count.data.1 <- data.frame(res.low.ref.with.counts,ensembl_gene_id = rownames(res.low.ref.with.counts))
  count.data.2 <- merge(count.data.1,G_list,by="ensembl_gene_id",sort = FALSE)
  count.data.2[which(count.data.2$mgi_symbol == ""),]$mgi_symbol <- paste0(count.data.2[which(count.data.2$mgi_symbol == ""),]$ensembl_gene_id,"-No_match")
  data.byGSym.order <- count.data.2[with(count.data.2,order(padj,log2FoldChange)),]
  write.table(data.byGSym.order,file = file.path(output.file.dir,paste0(fileName,".txt")),append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F)
}

fileName <- "LowVsRef"
getGeneSymbolOut(res.low.ref.with.counts,G_list,fileName)

fileName <- "HighVsRef"
getGeneSymbolOut(res.high.ref.with.counts,G_list,fileName)

fileName <- "HighVsLow"
getGeneSymbolOut(res.high.low.with.counts,G_list,fileName)

# 
# thresLogFc = 1
# thresP = 0.05
# sig.index <- with(data.byGSym.order,which(abs(log2FoldChange) > thresLogFc & padj < thresP))
# data.byGSym.order.sig <- data.byGSym.order[sig.index,]
# 
# EnsemToGene <- data.frame(Geneid=data.byGSym.order.sig[,1],hgnc_symbol=data.byGSym.order.sig[,16])
# 
# colnames(EnsemToGene) <- c("Geneid","hgnc_symbol")
# 
# rawdata.tli.after.filter.by.cpm <- data.byGSym.order.sig[,8:11]
# rownames(rawdata.tli.after.filter.by.cpm) <- data.byGSym.order.sig[,16]
# 
# tli.colData <- data.frame(condition=rep(c("ASA","Vector"),c(2,2)))
# rownames(tli.colData) <- colnames(data.byGSym.order.sig[,8:11])
# 
# rld <- DESeqDataSetFromMatrix(
#   countData=round(rawdata.tli.after.filter.by.cpm,0),
#   colData=tli.colData,design= ~condition)
# 
# library(gplots)
# library(RColorBrewer)
# 
# normalized.data <- data.byGSym.order.sig[,12:15]
# rownames(normalized.data) <- data.byGSym.order.sig[,16]
# 
# mydata <- heatmap_2group(normalized.data,data.byGSym.order.sig$padj,data.byGSym.order.sig$log2FoldChange,"ASAvsVector",output.file.dir,cutoffs=c(0.05,1),EnsemToGene=EnsemToGene,numCluster=2,rld=rld,which.genes=c("GCNT1","GCNT2"))
# 
# out.file.name <- "ASA-Vector"
# rowName <- rownames(mydata$mydata)
# geneSymbol <- EnsemToGene[match(rowName,EnsemToGene$hgnc_symbol),]$hgnc_symbol
# selected.gene <- c("GCNT1","GCNT2")
# groups <- c("ASA", "Vector")
# library(ComplexHeatmap)
# library(cowplot)
# heatmap.3(mydata$mydata,rowName = rowName,geneSymbol= geneSymbol,selected.gene = selected.gene,groups = groups,output.file.dir,out.file.name,g1Name="ASA",g2Name="Vector",g1N=2,g2N=2,w=5,h=8)

save.image(file = file.path(output.file.dir,"DE.RData"))
savehistory(file = file.path(output.file.dir,"DE.Rhistory"))
