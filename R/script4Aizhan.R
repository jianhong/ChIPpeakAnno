#!/usr/bin/env Rscript

output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output4Aizhan"

if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}

cm.unique.only.1.anno <- getAnnotatedGene(cm.unique.only.1,"Mm")

plotFPKM4AizhanNew <- function(){

  YYY.Aizhan <- getFPKM4DiffSet(cm.unique.only.1.anno,rna.seq.data.from.Aizhan)
  
  # /Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Aizhan_boxplot_2_5_2019/F121.only-XL.only-F121.and.XL-XL_NAD-NAD-nonNAD-wholeGenome.png
  YYY.Aizhan$SetName
  # 
  # /Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/ESCs_boxplot/F121.only-XL.only-F121.and.XL-NAD-nonNAD-wholeGenome.png
  # YYY.ESCs
  YYY.ESCs <- getFPKM4DiffSet(cm.unique.only.1.anno,rna.seq.data.ESCs)
  
  # /Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/ES1_boxplot/F121.only-XL.only-F121.and.XL-NAD-nonNAD-wholeGenome.png
  # 
  
  YYY.ES1 <- getFPKM4DiffSet(cm.unique.only.1.anno,rna.seq.data.ES1.rpkm)
  
  YYY.ES1
  # /Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/MEF_boxplot_2_6_2019/F121.only-XL.only-F121.and.XL-XL_NAD-NAD-nonNAD-wholeGenome.png
  # 
  YYY1 <- getFPKM4DiffSet(cm.unique.only.1.anno,fpkm.value)
  # /Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/MEF_boxplot_using_GSE90894_2_13_2019/F121.only-XL.only-F121.and.XL-XL_NAD-NAD-nonNAD-wholeGenome.png
  
  # YYY1.GSE90894.mef
  YYY1.GSE90894.mef <- getFPKM4DiffSet(cm.unique.only.1.anno,GSE90894.rnaSeq.data.2)
  # /Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/F121_9_RNA_Seq_data/GSM1418813_ES1_counts.txt
  # /Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/F121_9_RNA_seq_data/GSE90894_RPKM_mRNAseq_table.xlsx
  # /Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/F121_9_RNA_seq_data/GSM1621026_WT_fpkm.csv
  
  getFPKMTableAndFig4SetPeaks <- function(interest.peak) {
    
    # interest.peak  <- "XL_NAD"
    
    YYY1.GSE90894.mef.XL_NAD <- YYY1.GSE90894.mef[YYY1.GSE90894.mef$SetName==interest.peak,]
    
    YYY1.XL_NAD <- YYY1[YYY1$SetName==interest.peak,]
    
    YYY.ES1.XL_NAD <- YYY.ES1[YYY.ES1$SetName==interest.peak,]
    
    YYY.ESCs.XL_NAD <- YYY.ESCs[YYY.ESCs$SetName==interest.peak,]
    
    YYY.Aizhan.XL_NAD <- YYY.Aizhan[YYY.Aizhan$SetName==interest.peak,]
    
    gene.com <- intersect(YYY1.GSE90894.mef.XL_NAD$GeneName,YYY1.XL_NAD$GeneName)
    
    gene.com.1 <- intersect(gene.com,YYY.ES1.XL_NAD$GeneName)
    
    gene.com.2 <- intersect(gene.com.1,YYY.ESCs.XL_NAD$GeneName)
    
    gene.com.3 <- intersect(gene.com.2,YYY.Aizhan.XL_NAD$GeneName)
    
    YYY1.GSE90894.mef.XL_NAD.subset.gene <- YYY1.GSE90894.mef.XL_NAD[match(gene.com.3,YYY1.GSE90894.mef.XL_NAD$GeneName),]
    
    YYY1.XL_NAD.subset.gene <- YYY1.XL_NAD[match(gene.com.3,YYY1.XL_NAD$GeneName),]
    
    YYY.ES1.XL_NAD.subset.gene <- YYY.ES1.XL_NAD[match(gene.com.3,YYY.ES1.XL_NAD$GeneName),]
    
    YYY.ESCs.XL_NAD.subset.gene <- YYY.ESCs.XL_NAD[match(gene.com.3,YYY.ESCs.XL_NAD$GeneName),]
    
    YYY.Aizhan.XL_NAD.subset.gene <- YYY.Aizhan.XL_NAD[match(gene.com.3,YYY.Aizhan.XL_NAD$GeneName),]
    
    YYY1.GSE90894.mef.XL_NAD.subset.gene.1 <- cbind(RnaSeq=rep("GSE90894.mef",dim(YYY1.GSE90894.mef.XL_NAD.subset.gene)[1]),YYY1.GSE90894.mef.XL_NAD.subset.gene)
    
    YYY1.XL_NAD.subset.gene.1 <- cbind(RnaSeq=rep("Collas_GSM1621026.mef",dim(YYY1.XL_NAD.subset.gene)[1]),YYY1.XL_NAD.subset.gene)
    
    YYY.ES1.XL_NAD.subset.gene.1 <- cbind(RnaSeq=rep("Lowe_mESC",dim(YYY.ES1.XL_NAD.subset.gene)[1]),YYY.ES1.XL_NAD.subset.gene)
    
    YYY.ESCs.XL_NAD.subset.gene.subset.gene.1 <- cbind(RnaSeq=rep("GSE90894.ESCs",dim(YYY.ESCs.XL_NAD.subset.gene)[1]),YYY.ESCs.XL_NAD.subset.gene)
    
    YYY.Aizhan.XL_NAD.subset.gene.1 <- cbind(RnaSeq=rep("F121-9_RNAseq",dim(YYY.Aizhan.XL_NAD.subset.gene)[1]),YYY.Aizhan.XL_NAD.subset.gene)
    
    sum.5.sets <- rbind(YYY1.GSE90894.mef.XL_NAD.subset.gene.1,YYY1.XL_NAD.subset.gene.1,YYY.ES1.XL_NAD.subset.gene.1,YYY.ESCs.XL_NAD.subset.gene.subset.gene.1,YYY.Aizhan.XL_NAD.subset.gene.1)
    
    jpeg(filename = file.path(output.file.dir,paste0("FPKM_5_data_set_",interest.peak,".jpeg")),width = 1000, height = 400)
    boxplot(as.numeric(sum.5.sets$FPKM)~sum.5.sets$RnaSeq)
    dev.off()
    
    sum.5.sets.short <- data.frame(GSE90894.mef=sum.5.sets[sum.5.sets$RnaSeq=="GSE90894.mef",2:4],Collas_GSM1621026.mef=sum.5.sets[sum.5.sets$RnaSeq=="Collas_GSM1621026.mef",4], Lowe_mESC=sum.5.sets[sum.5.sets$RnaSeq=="Lowe_mESC",4],GSE90894.ESCs=sum.5.sets[sum.5.sets$RnaSeq=="GSE90894.ESCs",4],F121_9_RNAseq=sum.5.sets[sum.5.sets$RnaSeq=="F121-9_RNAseq",4])
    
    colnames(sum.5.sets.short)[c(1,2)] = c("Peak_set","Gene_Name")
    
   # output.file.name <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/FPKM_5_data_set_XL_NAD.txt"
    
    output.file.name <- file.path(output.file.dir,paste0("FPKM_5_data_set_",interest.peak,".txt"))
    
    write.table(sum.5.sets.short,file = output.file.name,append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T)
    
  }
  
  interest.peak  <- "XL_NAD"
  getFPKMTableAndFig4SetPeaks(interest.peak)
  
  interest.peak  <- "XL.only"
  getFPKMTableAndFig4SetPeaks(interest.peak)
  
}

save.image(file = file.path(output.file.dir,"FPKM.RData"))
savehistory(file = file.path(output.file.dir,"FPKM.Rhistory"))

