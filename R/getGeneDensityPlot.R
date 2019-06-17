# bed.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Results_10_4_2018_gene_density"

# bed.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Results_10_7_2018_gene_density"
# getGeneDensityPlot(bed.dir)

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

load("/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4Aizhan_F121_9_avesig_filtered_NADsVenn/F121-9_box_plot/F121-9_box_plot.RData")

bed.4.gene.density.1.anno <- getAnnotatedGene(bed.4.gene.density.1,"Mm")
bed.4.gene.density.1.anno.YYY <- getFPKM4DiffSet(bed.4.gene.density.1.anno,rna.seq.data.from.Aizhan)

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

bed.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4AizhanGeneDensity"

if(!dir.exists(bed.dir)){dir.create(bed.dir,recursive = TRUE)}
ggexport(multi.page, width = 1000, height = 480,filename = file.path(bed.dir,"GeneDensity.png"))

write.table(df[,1:4],file = file.path(bed.dir,paste0("Density",".txt")),
            append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T)



getGeneDensityPlot2 <- function(bed.dir) {
  
  if(!dir.exists(bed.dir)){dir.create(bed.dir,recursive = TRUE)}
  # whole genome in MB
  whole.genome.MB <- sum(width(gr.whole.genome))/1000000
  
  names(XL.nonXL.subset) <- c("XL_Uniq","XL_ciLAD","XL_LAD","nonXL_Uniq","nonXL_ciLAD","nonXL_LAD")
  
  null <- lapply(1:length(XL.nonXL.subset), function(u,XL.nonXL.subset){
    
    x <- names(XL.nonXL.subset)[u] 
    export(XL.nonXL.subset[[u]],file.path(bed.dir,paste0(x,".bed")),format="BED")
    
  },XL.nonXL.subset)
  
  xl.nad.bed <- file.path(bed.dir,paste0(c("XL_Uniq","XL_ciLAD","XL_LAD"),".bed"))
  cmd = paste("bedops -u",paste(xl.nad.bed,collapse = " "),">",file.path(bed.dir,"XLNAD.bed"))
  system(cmd)
  cmd = paste("bedtools sort -i",file.path(bed.dir,"XLNAD.bed"),">",file.path(bed.dir,"XLNAD_sorted.bed"))
  system(cmd) 
  cmd = paste("bedtools merge -i",file.path(bed.dir,"XLNAD_sorted.bed"),">",file.path(bed.dir,"XLNAD_merged.bed"))
  system(cmd) 
  cmd = paste("bedops -d",ref.bed,file.path(bed.dir,"XLNAD_merged.bed"),">",file.path(bed.dir,"XL_nonNAD.bed"))
  system(cmd)
  
  
  peaks=read.table(file.path(bed.dir,"XLNAD_merged.bed"))
  colnames(peaks)[1:3]= c("chr","start","end")
  XLNAD=toGRanges(peaks)
  
  peaks=read.table(file.path(bed.dir,"XL_nonNAD.bed"))
  colnames(peaks)[1:3]= c("chr","start","end")
  XL_nonNAD=toGRanges(peaks)
  
  
  
  nonxl.nad.bed <- file.path(bed.dir,paste0(c("nonXL_Uniq","nonXL_ciLAD","nonXL_LAD"),".bed"))
  cmd = paste("bedops -u",paste(nonxl.nad.bed,collapse = " "),">",file.path(bed.dir,"nonXLNAD.bed"))
  system(cmd)
  cmd = paste("bedtools sort -i",file.path(bed.dir,"nonXLNAD.bed"),">",file.path(bed.dir,"nonXLNAD_sorted.bed"))
  system(cmd) 
  cmd = paste("bedtools merge -i",file.path(bed.dir,"nonXLNAD_sorted.bed"),">",file.path(bed.dir,"nonXLNAD_merged.bed"))
  system(cmd) 
  cmd = paste("bedops -d",ref.bed,file.path(bed.dir,"nonXLNAD_merged.bed"),">",file.path(bed.dir,"nonXL_nonNAD.bed"))
  system(cmd)
  
  peaks=read.table(file.path(bed.dir,"nonXLNAD_merged.bed"))
  colnames(peaks)[1:3]= c("chr","start","end")
  nonXLNAD=toGRanges(peaks)
  
  peaks=read.table(file.path(bed.dir,"nonXL_nonNAD.bed"))
  colnames(peaks)[1:3]= c("chr","start","end")
  nonXL_nonNAD=toGRanges(peaks)
  
  XLNAD.MB <- sum(width(XLNAD))/1000000
  XL_nonNAD.MB <- sum(width(XL_nonNAD))/1000000
  
  nonXLNAD.MB <- sum(width(nonXLNAD))/1000000
  nonXL_nonNAD.MB <- sum(width(nonXL_nonNAD))/1000000
  
  
  # export(re.out[[6]],file.path(bed.dir,"nonXL.bed"),format="BED")
  # export(re.out[[8]],file.path(bed.dir,"XL.bed"),format="BED")
  # 
  # nad.bed <- paste(file.path(bed.dir,"nonXL.bed"),file.path(bed.dir,"XL.bed"),sep= " ")
  # cmd = paste("bedops -u",nad.bed,">",file.path(bed.dir,"NAD.bed"))
  # system(cmd)
  # cmd = paste("bedtools sort -i",file.path(bed.dir,"NAD.bed"),">",file.path(bed.dir,"NAD_sorted.bed"))
  # system(cmd) 
  # cmd = paste("bedtools merge -i",file.path(bed.dir,"NAD_sorted.bed"),">",file.path(bed.dir,"NAD_merged.bed"))
  # system(cmd) 
  # peaks=read.table(file.path(bed.dir,"NAD_merged.bed"))
  # colnames(peaks)[1:3]= c("chr","start","end")
  # NAD=toGRanges(peaks)
  # 
  # # NAD.before <- unique(c(re.out[[6]],re.out[[8]]))
  # 
  # XLNAD <- re.out[[8]]
  # nonXLNAD <- re.out[[6]]
  #   
  # # export(NAD,file.path(bed.dir,"NAD.bed"),format="BED")
  # 
  # 
  # cmd = paste("bedops -d",ref.bed,file.path(bed.dir,"NAD_merged.bed"),">",file.path(bed.dir,"nonNAD.bed"))
  # 
  # system(cmd)
  # 
  # peaks=read.table(file.path(bed.dir,"nonNAD.bed"))
  # 
  # colnames(peaks)[1:3]= c("chr","start","end")
  # nonNAD=toGRanges(peaks)
  # 
  # #mcols(NAD) <- NULL
  # #nonNAD <- setdiff(gr.whole.genome,NAD)
  # 
  # XLNAD.MB <- sum(width(XLNAD))/1000000
  # nonXLNAD.MB <- sum(width(nonXLNAD))/1000000
  # 
  # NAD.MB <- sum(width(NAD))/1000000
  # nonNAD.MB <- sum(width(nonNAD))/1000000
  
  mbl <- lapply(1:length(XL.nonXL.subset),function(u,XL.nonXL.subset){
    
    mb <- sum(width(XL.nonXL.subset[[u]]))/1000000
    mb
  },XL.nonXL.subset)
  names(mbl) <- names(XL.nonXL.subset)
  #names(mbl) <- unique(as.character(YYYY$SetName))[1:6]
  
  mbl2 <- do.call(rbind,mbl)
  colnames(mbl2) <- "MB"
  
  mbl3 <- as.data.frame(c(XLNAD.MB,XL_nonNAD.MB,nonXLNAD.MB,nonXL_nonNAD.MB,whole.genome.MB))
  
  
  
  row.names(mbl3) <- c("XLNAD.MB","XL_nonNAD.MB","nonXLNAD.MB","nonXL_nonNAD.MB","whole.genome.MB")
  colnames(mbl3) <- "MB"
  mbl4 <- rbind.data.frame(mbl2,mbl3)
  
  YYYYY <- YYYY
  
  YYYY.wholeGenome <- YYYY[which(YYYY$SetName %in% c("wholeGenome")),]
  
  YYYY.XLNAD <- YYYY[which(YYYY$SetName %in% c("XL_Uniq","XL_ciLAD","XL_LAD")),]
  YYYY.XLNAD.3.class <- YYYY.XLNAD
  YYYY.XLNAD$SetName <- "XLNAD" 
  
  YYYY.XL_nonNAD <- YYYY.wholeGenome[-which(YYYY.wholeGenome$GeneName %in% YYYY.XLNAD$GeneName),]
  YYYY.XL_nonNAD$SetName <- "XL_nonNAD" 
  
  YYYY.nonXLNAD <- YYYY[which(YYYY$SetName %in% c("nonXL_Uniq","nonXL_ciLAD","nonXL_LAD")),]
  YYYY.nonXLNAD.3.class <- YYYY.nonXLNAD
  YYYY.nonXLNAD$SetName <- "nonXLNAD" 
  
  YYYY.nonXL_nonNAD <- YYYY.wholeGenome[-which(YYYY.wholeGenome$GeneName %in% YYYY.nonXLNAD$GeneName),]
  YYYY.nonXL_nonNAD$SetName <- "nonXL_nonNAD"
  
  YYYY <- rbind.data.frame(YYYY.XLNAD.3.class,YYYY.XLNAD,YYYY.XL_nonNAD,YYYY.nonXLNAD.3.class,YYYY.nonXLNAD,YYYY.nonXL_nonNAD,YYYY.wholeGenome)
  
  YYYY$SetName <- droplevels(YYYY$SetName)
  
  png(file = file.path(bed.dir,paste0(paste(levels(YYYY$SetName),collapse = "-"),".png")),width = 1500, height = 480)
  boxplot(log10(as.numeric(as.character(YYYY$FPKM))+1)~YYYY$SetName,ylab = "log10(FPKM+1)")
  dev.off()
  
  
  num.gene <- lapply(unique(as.character(YYYY$SetName)), function(u,YYYY){
    num.gene <- length(unique(as.character(YYYY[which(YYYY$SetName==u),]$GeneName)))
    num.gene
  },YYYY)
  names(num.gene) <- unique(as.character(YYYY$SetName))
  num.gene.1 <- do.call(rbind,num.gene)
  colnames(num.gene.1) <- "num"
  
  row.names(mbl4) <- gsub(".MB","",row.names(mbl4))
  row.names(mbl4) <- gsub("whole.genome","wholeGenome",row.names(mbl4))
  
  mbl5 <- merge(mbl4,num.gene.1,by=0)
  
  #XLNAD.num.gene <- length(unique(as.character(YYYY[which(YYYY$SetName %in% c("XL_Uniq","XL_ciLAD","XL_LAD")),]$GeneName)))
  #nonXLNAD.num.gene <- length(unique(as.character(YYYY[which(YYYY$SetName %in% c("nonXL_Uniq","nonXL_ciLAD","nonXL_LAD")),]$GeneName)))
  
  #XLNAD.nonXLNAD.num.gene <- as.data.frame(c(XLNAD.num.gene,nonXLNAD.num.gene))
  #colnames(XLNAD.nonXLNAD.num.gene) <- "num"
  #row.names(XLNAD.nonXLNAD.num.gene) <- c("XLNAD","nonXLNAD")
  
  #num.gene.1 <- rbind(num.gene.1,XLNAD.nonXLNAD.num.gene)
  
  
  #mbl4[rownames(mbl4)=="XLNAD",] <- sum(mbl4[rownames(mbl4) %in% c("XL_Uniq","XL_ciLAD","XL_LAD"),])
  
  #mbl4[rownames(mbl4)=="nonXLNAD",] <- sum(mbl4[rownames(mbl4) %in% c("nonXL_Uniq","nonXL_ciLAD","nonXL_LAD"),])
  
  mb.gene.num <- cbind(mbl5$MB,mbl5$num,mbl5$num/mbl5$MB,mbl5$Row.names)
  
  #mb.gene.num <- data.frame(mb.gene.num,Category=row.names(mb.gene.num))
  
  colnames(mb.gene.num) <- c("MB","num.gene","Density","Category")
  
  DF <- as.data.frame(mb.gene.num)
  DF$MB <- as.numeric(as.character(DF$MB))
  DF$num.gene <- as.numeric(as.character(DF$num.gene))
  DF$Density <- as.numeric(as.character(DF$Density))
  
  #create a colour scheme based on grouping variable 'zone'
  Category<- as.character(DF$Category)
  color.codes <- rainbow(length(Category))
  
  c.1 <- color.codes[which(Category %in% c("nonXL_Uniq","XL_Uniq"))][1]
  
  color.codes[which(Category %in% c("nonXL_Uniq","XL_Uniq"))] <- c.1
  
  c.2 <- color.codes[which(Category %in% c("nonXL_ciLAD","XL_ciLAD"))][1]
  color.codes[which(Category %in% c("nonXL_ciLAD","XL_ciLAD"))] <- c.2
  
  c.3 <- color.codes[which(Category %in% c("nonXL_LAD" ,"XL_LAD"))][1]
  color.codes[which(Category %in% c("nonXL_LAD" ,"XL_LAD"))] <- c.3
  
  c.4 <- color.codes[which(Category %in% c("nonXLNAD" ,"XLNAD"))][1]
  color.codes[which(Category %in% c("nonXLNAD" ,"XLNAD"))] <- c.4
  
  #color.codes<-as.character(c("#3399FF", "#FF0000","#9633FF","#3399FF", "#FF0000","#9633FF","#0000FF","#C71585","#FFC0CB","#FF0000FF","#00FFFFFF"))
  
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
  
}
