#!/usr/bin/env Rscript

# 1) F121-9 NADs
# Dropbox/pub_old_galaxy/F121-9/avesig filtered NADs/F121_9_all_1.5FC_F200_minp_adjpVal_avesig_greater_than_1.7_NAD_greater_50kb.bed

F121.9.NADs <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/avesig filtered NADs/F121_9_all_1.5FC_F200_minp_adjpVal_avesig_greater_than_1.7_NAD_greater_50kb.bed"

# 2) F121-9 Type I NADs
# Dropbox\pub_old_galaxy\F121-9\overlap analysis\F121-9 vs CLAD vs ciLAD\common-in-F121_9_and_cLAD.bed

F121.9.Type.I.NADs <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/overlap analysis/F121-9 vs CLAD vs ciLAD/common-in-F121_9_and_cLAD.bed"

# 3) F121-9 Type II NADs
# Dropbox\pub_old_galaxy\F121-9\overlap analysis\F121-9 vs CLAD vs ciLAD\common-in-F121_9_and_ciLAD.bed

F121.9.Type.II.NADs <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/overlap analysis/F121-9 vs CLAD vs ciLAD/common-in-F121_9_and_ciLAD.bed"

# 4) Shared NADs
# Dropbox\pub_old_galaxy\F121-9\overlap analysis\F121-9 vs MEF XL NADs\common-in-F121_9_and-XL_MEF.bed

shared.NADs <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/overlap analysis/F121-9 vs MEF XL NADs/common-in-F121_9_and-XL_MEF.bed"

# 5) F121-9 specific NADs
# Dropbox\pub_old_galaxy\F121-9\overlap analysis\F121-9 vs MEF XL NADs\unique-F121_9.bed

F121.9.specific.NADs <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/overlap analysis/F121-9 vs MEF XL NADs/unique-F121_9.bed"

# 6) MEF specific NADs
# Dropbox\pub_old_galaxy\F121-9\overlap analysis\F121-9 vs MEF XL NADs\unique-XL_MEF.bed

MEF.specific.NADs <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/overlap analysis/F121-9 vs MEF XL NADs/unique-XL_MEF.bed"

# 7) cLAD
# Dropbox/NADfinder/BedFiles/cLAD.mm10 bed with header.bed

cLAD <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/BedFiles/cLAD.mm10 bed with header.bed"

# 8) ciLAD
# Dropbox/NADfinder/BedFiles/ciLAD.mm10 bed with header.bed

ciLAD <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/BedFiles/ciLAD.mm10 bed with header.bed"

# 9) Early replication timing: all regions above zero in this bedgraph 
# Dropbox (UMass Medical School)\pub_old_galaxy\F121-9\published data\GSE95091_F121-9_nesc_r1_R_T_.bedGraph

Early.replication.timing <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/published data/GSE95091_F121-9_nesc_r1_R_T_.bedGraph"

# 10) Late replication timing: all regions below zero in this bedgraph
# Dropbox (UMass Medical School)\pub_old_galaxy\F121-9\published data\GSE95091_F121-9_nesc_r1_R_T_.bedGraph

Late.replication.timing <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/published data/GSE95091_F121-9_nesc_r1_R_T_.bedGraph"

H3K27me3 <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/H3K27me3/GSM2416833_WT_2i_ChIP_H3K27me3-domains-1.bed"

F121_9_nonNADs <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4_F121_9_Type_I_NADs_F121_9_Type_II_NADs_F121_9_NADs_F121_9_nonNADs_box_plot/F121_9_nonNADs.bed"


output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4AizhanJaccardAnalysis"

bed.files <- c(F121.9.NADs,F121.9.Type.I.NADs,F121.9.Type.II.NADs,shared.NADs,F121.9.specific.NADs,MEF.specific.NADs,cLAD,ciLAD,Early.replication.timing,Late.replication.timing,H3K27me3,F121_9_nonNADs)

bed.in<-lapply(1:length(bed.files),function(u,bed.files){
  
  x <- bed.files[[u]]
  
  if(u == 2|u==3|u==4|u==5|u==6|u==9|u==10|u==11|u==12){peaks=read.table(x,header=F)}else{peaks=read.table(x,header=F,skip=1)}
  
  colnames(peaks)[1:3]= c("chr","start","end")
  peaks=toGRanges(peaks)
  peaks
},bed.files)

names(bed.in) <- c("F121.9.NADs","F121.9.Type.I.NADs","F121.9.Type.II.NADs","shared.NADs","F121.9.specific.NADs","MEF.specific.NADs","cLAD","ciLAD","Early.replication.timing","Late.replication.timing","H3K27me3","F121_9_nonNADs")

Ert <- bed.in[[9]][which(bed.in[[9]]$V4>0)]
Lrt <- bed.in[[10]][which(bed.in[[10]]$V4<0)]

bed.in[[9]] <- Ert
bed.in[[10]] <- Lrt


names(bed.in[c(12,3,1,2)])

bed.4.gene.density <- bed.in[c(12,3,1,2)]
gr.whole.genome <- GRanges(seqinfo(BSgenome.Mmusculus.UCSC.mm10))

levels.to.keep <- levels(seqnames(gr.whole.genome))[1:21]
w.gr.simple<- keepSeqlevels(gr.whole.genome,levels.to.keep,pruning.mode="coarse")

genome(seqinfo(w.gr.simple)) <- "GRCm38"
bed.4.gene.density.1 <- c(bed.4.gene.density[1],gr.whole.genome=w.gr.simple,bed.4.gene.density[c(2:4)])

Gr2DfWithMB <- function(XL.nonXL.subset) {
  mbl <- lapply(1:length(XL.nonXL.subset),function(u,XL.nonXL.subset){
    mb <- sum(width(XL.nonXL.subset[[u]]))/1000000
    mb
  },XL.nonXL.subset)
  names(mbl) <- names(XL.nonXL.subset)
  mbl2 <- do.call(rbind,mbl)
  colnames(mbl2) <- "MB"
  mbl2 
}

bed.4.gene.density.1.df <- Gr2DfWithMB(bed.4.gene.density.1)

bed.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4AizhanGeneDensity"

rna.seq.data.from.Aizhan <- readRDS(file=file.path(bed.dir,"rnaSeq_data_for_F121_9_FPKM.rds"))

#load("/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4Aizhan_F121_9_avesig_filtered_NADsVenn/F121-9_box_plot/F121-9_box_plot.RData")

bed.4.gene.density.1.anno <- getAnnotatedGene(bed.4.gene.density.1,"Mm")
bed.4.gene.density.1.anno.YYY <- getFPKM4DiffSet(bed.4.gene.density.1.anno,rna.seq.data.from.Aizhan)

getFPKMboxplot <- function(bed.dir,YYYY) {
  
  png(file = file.path(bed.dir,paste0(paste(levels(YYYY$SetName),collapse = "-"),".png")),width = 1500, height = 480)
  boxplot(log10(as.numeric(as.character(YYYY$FPKM))+1)~YYYY$SetName,ylab = "log10(FPKM+1)")
  dev.off()

}

getFPKMboxplot(bed.dir,bed.4.gene.density.1.anno.YYY)

getNumOfGene4Sets <- function(YYYY) {
  num.gene <- lapply(unique(as.character(YYYY$SetName)), function(u,YYYY){
    num.gene <- length(unique(as.character(YYYY[which(YYYY$SetName==u),]$GeneName)))
    num.gene
  },YYYY)
  names(num.gene) <- unique(as.character(YYYY$SetName))
  num.gene.1 <- do.call(rbind,num.gene)
  colnames(num.gene.1) <- "num"
  num.gene.1
}

bed.4.gene.density.1.anno.YYY.gene.num <- getNumOfGene4Sets(bed.4.gene.density.1.anno.YYY)

mbl5 <- merge(bed.4.gene.density.1.df,bed.4.gene.density.1.anno.YYY.gene.num,by=0)

mb.gene.num <- cbind(mbl5$MB,mbl5$num,mbl5$num/mbl5$MB,mbl5$Row.names)

colnames(mb.gene.num) <- c("MB","num.gene","Density","Category")

DF <- as.data.frame(mb.gene.num)
DF$MB <- as.numeric(as.character(DF$MB))
DF$num.gene <- as.numeric(as.character(DF$num.gene))
DF$Density <- as.numeric(as.character(DF$Density))

#create a colour scheme based on grouping variable 'zone'
Category<- as.character(DF$Category)
color.codes <- rainbow(length(Category))

df2=data.frame(Category, color.codes)

# merge color specifications with data
df <-merge(DF,df2, by=("Category"), all.x=TRUE, all.y=TRUE)

df$Category <- factor(df$Category, levels = df$Category[order(df$Density,decreasing = T)])

sp <- ggplot(df, aes(x = Category, y = Density, fill = color.codes)) +
  geom_bar(stat="identity",position = "dodge")  +
  theme(plot.title=element_text(hjust=0.5)) +  
  ylab("Gene density (# of genes/MB)")  +  theme(legend.position="none")

spp <- list(sp=sp)
multi.page <- ggarrange(plotlist=spp,nrow = 1, ncol = 1)

if(!dir.exists(bed.dir)){dir.create(bed.dir,recursive = TRUE)}
ggexport(multi.page, width = 1000, height = 480,filename = file.path(bed.dir,"GeneDensity.png"))

write.table(df[,1:4],file = file.path(bed.dir,paste0("Density",".txt")),
            append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T)

saveRDS(bed.in,file = file.path(bed.dir,paste0(paste(names(bed.in),collapse = "-"),".rds")))
save.image(file = file.path(bed.dir,paste0(basename(bed.dir),".RData")))
savehistory(file = file.path(bed.dir,paste0(basename(bed.dir),".Rhistory")))

