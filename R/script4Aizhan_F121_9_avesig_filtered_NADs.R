#!/usr/bin/env Rscript

input.file <- "~/Aimin/DropboxUmass/NADfinder/Aimin/BamFile_From_Aizhan/count2.txt"
input.file.dir <- "~/Aimin/DropboxUmass/pub_old_galaxy/F121-9/avesig filtered NADs"
output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output4Aizhan_F121_9_avesig_filtered_NADs"
if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}

getBedFiles5 <- function(input.file.dir) {
  file.name.bedgraph <- list.files(input.file.dir,pattern="*.bed$",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
  
  bed.in<-lapply(file.name.bedgraph,function(u){
    
    peaks=read.table(u,skip=1,header = F)
    colnames(peaks)[1:3]= c("chr","start","end")
    peaks=toGRanges(peaks)
    peaks
  })
  
  names(bed.in) <- gsub(" ","_",tools::file_path_sans_ext(file.name.bedgraph))
  names(bed.in) <- gsub("_mm10_copy","",names(bed.in))
  bed.in
  
}

input.bed <- getBedFiles5(input.file.dir)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("EnsDb.Mmusculus.v75", version = "3.8")
library(EnsDb.Mmusculus.v75)

cm.unique.only.1.anno <- getAnnotatedGene(input.bed,"Mm")

count2dataTable <- function(input.file) {
  
  count.data<- read.table(input.file,header= T)
  mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
  genes <- count.data$Geneid
  
  count.data.1 <- data.frame(count.data,ensembl_gene_id = genes)
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol"),values=count.data.1$ensembl_gene_id,mart= mart)
  count.data.2 <- merge(count.data.1,G_list,by="ensembl_gene_id")
  count.data.2[which(count.data.2$mgi_symbol %in% ""),]$mgi_symbol <- paste0(count.data.2[which(count.data.2$mgi_symbol %in% ""),]$ensembl_gene_id,"-No_match")
  
  countToFpkm <- function(counts, effLen)
  {
    N <- sum(counts)
    exp( log(counts) + log(1e9) - log(effLen) - log(N) )
  }
  
  fpkm <- data.frame(fpkm.8=countToFpkm(count.data.2[,8],count.data.2[,7]),fpkm.9=countToFpkm(count.data.2[,9],count.data.2[,7]))
  
  rna.seq.data.from.Aizhan <- data.frame(gene.name=count.data.2$mgi_symbol,fpkm=apply(fpkm,1,mean))
  rna.seq.data.from.Aizhan
}

rna.seq.data.from.Aizhan <- count2dataTable(input.file)

plotFPKM4AizhanNew <- function(cm.unique.only.1.anno,rna.seq.data.from.Aizhan){

  YYY.Aizhan <- getFPKM4DiffSet(cm.unique.only.1.anno,rna.seq.data.from.Aizhan)
  YYY.Aizhan$SetName <- gsub("<","_less_",YYY.Aizhan$SetName)
  YYY.Aizhan$SetName <- gsub(">","_greater_",YYY.Aizhan$SetName)
  YYY.Aizhan$SetName <- basename(YYY.Aizhan$SetName)
  YYY.Aizhan$SetName <- gsub("F121_9_all_1.5FC_F200_minp_adjPval_","",YYY.Aizhan$SetName)
  
  YYY.Aizhan.rm.outlier <- YYY.Aizhan[-which(as.numeric(as.character(YYY.Aizhan$FPKM))>=200),]
  
  png(file = file.path(output.file.dir,paste0(paste(unique(YYY.Aizhan$SetName),collapse = "-"),"-FPKM-",".png")),width = 1500, height = 480)
  boxplot(as.numeric(as.character(YYY.Aizhan.rm.outlier$FPKM))~YYY.Aizhan.rm.outlier$SetName,ylab = "FPKM")
  dev.off()
  
  png(file = file.path(output.file.dir,paste0(paste(unique(YYY.Aizhan$SetName),collapse = "-"),".png")),width = 1500, height = 480)
  boxplot(log10(as.numeric(as.character(YYY.Aizhan$FPKM))+1)~YYY.Aizhan$SetName,ylab = "log10(FPKM+1)")
  dev.off()
  
  output.file.name <- file.path(output.file.dir,paste0("FPKM_4_data_set_",".txt"))
  write.table(YYY.Aizhan,file = output.file.name,append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T)

}

plotFPKM4AizhanNew(cm.unique.only.1.anno,rna.seq.data.from.Aizhan)

save.image(file = file.path(output.file.dir,"FPKM.RData"))
savehistory(file = file.path(output.file.dir,"FPKM.Rhistory"))