#!/usr/bin/env Rscript


# We are now ready to start making what we hope will be the final set of figures for the paper submission.
# First, can you please do GO analysis for the following NADs:
#   
#   1)      F121-9 specific NADs:
#   Dropbox\pub_old_galaxy\F121-9\overlap analysis\F121-9 vs MEF XL NADs\unique-F121_9.bed
# 

unique.F121_9 <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/overlap analysis/F121-9 vs MEF XL NADs/unique-F121_9.bed"

# 2)      MEF specific NADs:
#   Dropbox\pub_old_galaxy\F121-9\overlap analysis\F121-9 vs MEF XL NADs\unique-XL_MEF.bed

unique.XL_MEF <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/overlap analysis/F121-9 vs MEF XL NADs/unique-XL_MEF.bed"

# 3)      Shared NADs:
#   Dropbox\pub_old_galaxy\F121-9\overlap analysis\F121-9 vs MEF XL NADs\common-in-F121_9_and-XL_MEF.bed
#

common.in.F121_9_and_XL_MEF <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/overlap analysis/F121-9 vs MEF XL NADs/common-in-F121_9_and-XL_MEF.bed"

# 4)      cLAD and F121-9 overlap:
#   Dropbox\pub_old_galaxy\F121-9\overlap analysis\F121-9 vs CLAD vs ciLAD\common-in-F121_9_and_cLAD.bed
#

common.in.F121_9_and_cLAD <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/overlap analysis/F121-9 vs CLAD vs ciLAD/common-in-F121_9_and_cLAD.bed"

# 5)      ciLAD and F121-9 overlap:
#   Dropbox\pub_old_galaxy\F121-9\overlap analysis\F121-9 vs CLAD vs ciLAD\common-in-F121_9_and_ciLAD.bed
#

common.in.F121_9_and_ciLAD <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/overlap analysis/F121-9 vs CLAD vs ciLAD/common-in-F121_9_and_ciLAD.bed"

input.files.new <- c(unique.F121_9,unique.XL_MEF,common.in.F121_9_and_XL_MEF,common.in.F121_9_and_cLAD,common.in.F121_9_and_ciLAD)

bed.in.new<-lapply(1:length(input.files.new),function(u,input.files.new){
  
  x <- input.files.new[[u]]
  
  if(u <= 5){peaks=read.table(x,header=F)}else{peaks=read.table(x,header=F,skip=1)}
  
  colnames(peaks)[1:3]= c("chr","start","end")
  peaks=toGRanges(peaks)
  peaks
},input.files.new)

names(bed.in.new) <- basename(input.files.new)
names(bed.in.new) <- c("unique_F121_9","unique_XL_MEF","F121_9_and_XL_MEF","F121_9_and_cLAD","F121_9_and_ciLAD")
bed.in.new.anno <- getAnnotatedGene(bed.in.new,"Mm")

output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/overlap analysis/Output4GOAnalysis"
genome <- "Mm"
bed.in.new.1.GO <- peaksToEnrichedGO(bed.in.new,output.file.dir,genome)

library(GO.db)
xx <- as.list(GOTERM)
foo <- function(x) c(GOID(x), Term(x),Definition(x), Ontology(x))
gomat <- t(sapply(xx, foo, simplify=TRUE))
gomat <- as.data.frame(gomat)
colnames(gomat) <- c("ID","Term","Definition","Ontology")
batchDraw4GO(output.file.dir,gomat)

save.image(file = file.path(output.file.dir,"GOAnalysis.RData"))
savehistory(file = file.path(output.file.dir,"GOAnalysis.Rhistory"))
