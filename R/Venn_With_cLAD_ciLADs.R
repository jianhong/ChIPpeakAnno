
output.file.dir <- "/Users/aiminyan/Aimin/DropboxUMass/NADfinder/Aimin/Venn_With_cLAD_ciLADs"

H3K27me3 <- "/Users/aiminyan/Aimin/DropboxUMass/NADfinder/BedFiles/GSM1621022_WT_H3K27me3_peaks.bed"

cLAD <- "/Users/aiminyan/Aimin/DropboxUMass/pub_old_galaxy/Kaufman_Lab_shared_files/cLAD.mm10 bed with header.bed"
  
XL_MEF <- "/Users/aiminyan/Aimin/DropboxUMass/NADfinder/BedFiles/Aizhan F121-9 comparisons/XL_MEF_1.5FC_minp_countF200_auto_q_less_than_e-3_chrX_q_less_than_e-2_adjPval copy - Copy.bed"

ciLAD <- "/Users/aiminyan/Aimin/DropboxUMass/pub_old_galaxy/Kaufman_Lab_shared_files/ciLAD.mm10 bed with header.bed"

bed.files <- c(H3K27me3,cLAD,XL_MEF,ciLAD)

bed.in<-lapply(1:length(bed.files),function(u,bed.files){
  
  x <- bed.files[[u]]
  
  if(u == 1){peaks=read.table(x,header=F)}else{peaks=read.table(x,header=F,skip=1)}
  
  colnames(peaks)[1:3]= c("chr","start","end")
  peaks=toGRanges(peaks)
  peaks
},bed.files)

names(bed.in) <- c("H3K27me3","cLAD","XL_MEF","ciLAD")

peak.index <- c(1,2,3)
name <- c("H3K27me3","cLAD","XL_MEF")
getCount4Venn(bed.in,peak.index,name,output.file.dir)

peak.index <- c(1,3,4)
name <- c("H3K27me3","XL_MEF","ciLAD")
getCount4Venn(bed.in,peak.index,name,output.file.dir)


