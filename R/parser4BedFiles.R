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

output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4AizhanJaccardAnalysis"

bed.files <- c(F121.9.NADs,F121.9.Type.I.NADs,F121.9.Type.II.NADs,shared.NADs,F121.9.specific.NADs,MEF.specific.NADs,cLAD,ciLAD,Early.replication.timing,Late.replication.timing)

bed.in<-lapply(1:length(bed.files),function(u,bed.files){
  
  x <- bed.files[[u]]
  
  if(u == 2|u==3|u==4|u==5|u==6|u==9|u==10){peaks=read.table(x,header=F)}else{peaks=read.table(x,header=F,skip=1)}
  
  colnames(peaks)[1:3]= c("chr","start","end")
  peaks=toGRanges(peaks)
  peaks
},bed.files)

names(bed.in) <- c("F121.9.NADs","F121.9.Type.I.NADs","F121.9.Type.II.NADs","shared.NADs","F121.9.specific.NADs","MEF.specific.NADs","cLAD","ciLAD","Early.replication.timing","Late.replication.timing")

Ert <- bed.in[[9]][which(bed.in[[9]]$V4>0)]
Lrt <- bed.in[[10]][which(bed.in[[10]]$V4<0)]

bed.in[[9]] <- Ert
bed.in[[10]] <- Lrt

if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}

getJC4aSetOfPeaks1 <- function(re.out.rt.mef.5.sets, output.file.dir) {
  d <- sapply(names(re.out.rt.mef.5.sets), function(x) sapply(names(re.out.rt.mef.5.sets), function(y) calculateJaccardCoefficient42Peaks(peaks1=re.out.rt.mef.5.sets[[y]], peaks2=re.out.rt.mef.5.sets[[x]])))
#  d <- d/100
  
  file.name <- paste0(length(re.out.rt.mef.5.sets),"_bed_")
  
  write.table(d,file = file.path(output.file.dir,paste0(file.name,"JC_similarity",".txt")),append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names=NA)
  d
}

d <- getJC4aSetOfPeaks1(bed.in, output.file.dir)

source("R/heatmap.3.R")

d.u <- unique(c(d))
d.u.order <- d.u[order(d.u)]
colors <- c(d.u.order)
my_palette <- redgreen(36)

file.name <- paste(names(bed.in),collapse = "_")
tiff(file.path(output.file.dir,paste0(file.name,"_heatmap4JC.tiff")))
heatmap.3(d,col=my_palette,Rowv = T,Colv=T,dendrogram="both",
          breaks=colors, density.info="none", trace="none", key=T,symm=T,symkey=F,symbreaks=T,margin=c(12, 13),key.xtickfun=function() {
            breaks <- parent.frame()$breaks
            return(list(
              at=parent.frame()$scale01(c(breaks[1],breaks[length(breaks)])),
              labels=c(as.character(breaks[1]),as.character(breaks[length(breaks)]))
            ))
          })
dev.off()

save.image(file = file.path(output.file.dir,paste0(basename(output.file.dir),".RData")))
savehistory(file = file.path(output.file.dir,paste0(basename(output.file.dir),".Rhistory")))