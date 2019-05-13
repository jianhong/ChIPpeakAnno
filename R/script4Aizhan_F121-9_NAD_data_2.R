#!/usr/bin/env Rscript

1) Here's the path to Avesig>1.7 file (I also removed NADs <50 kb):

Dropbox/pub_old_galaxy/F121-9/avesig filtered NADs/F121_9_all_1.5FC_F200_minp_adjpVal_avesig>1.7 NAD>50kb.bed

input.file1 <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/avesig filtered NADs/F121_9_all_1.5FC_F200_minp_adjpVal_avesig_greater_than_1.7_NAD_greater_50kb.bed"


2) ciLAD and cLAD are at:

Dropbox/NADfinder/BedFiles/cLAD.mm10 bed with header.bed
input.file2 <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/BedFiles/cLAD.mm10 bed with header.bed"

Dropbox/NADfinder/BedFiles/ciLAD.mm10 bed with header.bed

input.file3 <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/BedFiles/ciLAD.mm10 bed with header.bed"


3) MEF XL NADs:

Dropbox/NADfinder/BedFiles/Aizhan F121-9 comparisons/XL_MEF_1.5FC_minp_countF200_auto_q<e-3_chrX_q<e-2_adjPval copy.bed

input.file4 <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/BedFiles/Aizhan F121-9 comparisons/XL_MEF_1.5FC_minp_countF200_auto_q_less_than_e-3_chrX_q_less_than_e-2_adjPval copy - Copy.bed"


We also need to ask Aimin to have common or unique NADs to be >50 kb. 
Thanks,
Aizhan

input.files <- c(input.file1,input.file2,input.file3,input.file4)
bed.in<-lapply(input.files,function(u){
peaks=read.table(u,header=F,skip=1)
colnames(peaks)[1:3]= c("chr","start","end")
peaks=toGRanges(peaks)
peaks
})

names(bed.in) <- basename(input.files)

output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4Aizhan_F121_9_avesig_filtered_NADsVenn"

peak.index <- c(1,2,3,4)
name <- c("F121_9","cLAD","ciLAD","XL_MEF")
getCount4Venn(bed.in,peak.index,name,output.file.dir)

peak.index <- c(1,4)
name <- c("F121_9","XL_MEF")
getCount4Venn(bed.in,peak.index,name,output.file.dir)

peak.index <- c(1,2,3)
name <- c("F121_9","cLAD","ciLAD")
getCount4Venn(bed.in,peak.index,name,output.file.dir)

