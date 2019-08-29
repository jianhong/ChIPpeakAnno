#1) nonNAD
#Dropbox\NADfinder\Aimin\Output4_F121_9_Type_I_NADs_F121_9_Type_II_NADs_F121_9_NADs_F121_9_nonNADs_box_plot\F121_9_nonNADs.bed
F121_9.nonNAD <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4_F121_9_Type_I_NADs_F121_9_Type_II_NADs_F121_9_NADs_F121_9_nonNADs_box_plot/F121_9_nonNADs.bed"

#2) whole genome
#Dropbox\NADfinder\Aimin\Output4Aizhan_F121_9_avesig_filtered_NADsVenn\F121-9_box_plot\rnaSeq_data_for_F121_9_FPKM
whole.genome <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4Aizhan_F121_9_avesig_filtered_NADsVenn/F121-9_box_plot/rnaSeq_data_for_F121_9_FPKM.txt"

#3) Type II NAD
#Dropbox\pub_old_galaxy\F121-9\overlap analysis\unfiltered_only_F121_input_greater_50_kb\common-in-F121_9_all_1.5FC_F200_minp_adjpVal_avesig_greater_than_1.7_NAD_greater_50kb_bed-and-ciLAD_mm10 bed with header.bed.bed
Type.II.NAD <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/overlap analysis/unfiltered_only_F121_input_greater_50_kb/common-in-F121_9_all_1.5FC_F200_minp_adjpVal_avesig_greater_than_1.7_NAD_greater_50kb_bed-and-ciLAD_mm10 bed with header.bed.bed"

#4)F121-9 NAD
#Dropbox/pub_old_galaxy/F121-9/avesig filtered NADs/F121_9_all_1.5FC_F200_minp_adjpVal_avesig_greater_than_1.7_NAD_greater_50kb.bed
F121_9.NAD <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/avesig filtered NADs/F121_9_all_1.5FC_F200_minp_adjpVal_avesig_greater_than_1.7_NAD_greater_50kb.bed"

#5) Type I NAD
#Dropbox\pub_old_galaxy\F121-9\overlap analysis\unfiltered_only_F121_input_greater_50_kb\common-in-F121_9_all_1.5FC_F200_minp_adjpVal_avesig_greater_than_1.7_NAD_greater_50kb.bed-and-cLAD.mm10 bed with header.bed.bed
Type.I.NAD <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/overlap analysis/unfiltered_only_F121_input_greater_50_kb/common-in-F121_9_all_1.5FC_F200_minp_adjpVal_avesig_greater_than_1.7_NAD_greater_50kb.bed-and-cLAD.mm10 bed with header.bed.bed"

bed.files <- c(F121_9.nonNAD,Type.II.NAD,F121_9.NAD,Type.I.NAD)

bed.in<-lapply(1:length(bed.files),function(u,bed.files){
  
  x <- bed.files[[u]]
  
  if(u == 1|u==2|u==4){peaks=read.table(x,header=F)}else{peaks=read.table(x,header=F,skip=1)}
  
  colnames(peaks)[1:3]= c("chr","start","end")
  peaks=toGRanges(peaks)
  peaks
},bed.files)

names(bed.in) <- c("F121_9.nonNAD","Type.II.NAD","F121_9.NAD","Type.I.NAD")

bed.4.gene.density <- bed.in[c(1,2,3,4)]
gr.whole.genome <- GRanges(seqinfo(BSgenome.Mmusculus.UCSC.mm10))

levels.to.keep <- levels(seqnames(gr.whole.genome))[1:21]
w.gr.simple<- keepSeqlevels(gr.whole.genome,levels.to.keep,pruning.mode="coarse")

genome(seqinfo(w.gr.simple)) <- "GRCm38"
bed.4.gene.density.1 <- c(bed.4.gene.density,gr.whole.genome=w.gr.simple)

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

bed.4.gene.density.1.anno <- getAnnotatedGene(bed.4.gene.density.1,"Mm")

fileName <- names(bed.4.gene.density.1.anno)
  
library(GenomicRanges)
library(data.table)

null <- lapply(1:length(bed.4.gene.density.1.anno),function(u,bed.4.gene.density.1.anno,bed.dir){
  
  x <- as.data.table(bed.4.gene.density.1.anno[[u]])
  y <- names(bed.4.gene.density.1.anno)[u]
  con <- as.character(file.path(bed.dir,paste0(y,".csv")))
  print(con)
  write.table(x,con,quote = F,sep=",",row.names = F)
  
},bed.4.gene.density.1.anno,bed.dir)

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

# Way1 to save
ggsave(file.path(bed.dir,"GeneDensity_8_20_2019_dpi_300.png"),plot=sp,dpi=300)
       
# Way2 to save
library(ggpubr)
spp <- list(sp=sp)
multi.page <- ggarrange(plotlist=spp,nrow = 1, ncol = 1)

if(!dir.exists(bed.dir)){dir.create(bed.dir,recursive = TRUE)}
ggexport(multi.page, width = 1000, height = 480,filename = file.path(bed.dir,"GeneDensity_8_20_2019.png"))

save.image(file = file.path(bed.dir,"GeneListAnnotated2Bed.RData"))
savehistory(file = file.path(bed.dir,"GeneListAnnotated2Bed.Rhistory"))
