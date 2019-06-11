#!/usr/bin/env Rscript

# 1) Here's the path to Avesig>1.7 file (I also removed NADs <50 kb):
# 
# Dropbox/pub_old_galaxy/F121-9/avesig filtered NADs/F121_9_all_1.5FC_F200_minp_adjpVal_avesig>1.7 NAD>50kb.bed
# 
# input.file1 <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/avesig filtered NADs/F121_9_all_1.5FC_F200_minp_adjpVal_avesig_greater_than_1.7_NAD_greater_50kb.bed"
# 
# 
# 2) ciLAD and cLAD are at:
# 
# Dropbox/NADfinder/BedFiles/cLAD.mm10 bed with header.bed
# input.file2 <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/BedFiles/cLAD.mm10 bed with header.bed"
# 
# Dropbox/NADfinder/BedFiles/ciLAD.mm10 bed with header.bed
# 
# input.file3 <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/BedFiles/ciLAD.mm10 bed with header.bed"
# 
# 
# 3) MEF XL NADs:
# 
# Dropbox/NADfinder/BedFiles/Aizhan F121-9 comparisons/XL_MEF_1.5FC_minp_countF200_auto_q<e-3_chrX_q<e-2_adjPval copy.bed
# 
# input.file4 <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/BedFiles/Aizhan F121-9 comparisons/XL_MEF_1.5FC_minp_countF200_auto_q_less_than_e-3_chrX_q_less_than_e-2_adjPval copy - Copy.bed"
# 
# 
# We also need to ask Aimin to have common or unique NADs to be >50 kb. 
# Thanks,
# Aizhan

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
nonXL.3.subsets <- getUniquePeaks(ol.ciLAD.non_XL.MEF_LAD)
sum(width(nonXL.3.subsets[[3]]))

# 2) Overlap between F121-9 NADs and cLAD
ol.ciLAD.non_XL.MEF_LAD <- findOverlapsOfPeaks(bed.in[c(1,2)])
Aizhan.peaks.in.unique <- getPeaksInUnique(ol.ciLAD.non_XL.MEF_LAD)
Aizhan.peaks.in.overlap <- getPeaksInOverlap(ol.ciLAD.non_XL.MEF_LAD,1)
Aizhan.3.sets.of.peaks <- c(Aizhan.peaks.in.unique,Aizhan.peaks.in.overlap[[3]])
names(Aizhan.3.sets.of.peaks)[3] <- names(Aizhan.peaks.in.overlap)[3]
output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4Aizhan_F121_9_avesig_filtered_NADsVenn"
outGrl(Aizhan.3.sets.of.peaks,speci="Mm",output.file.dir)

# 3) Overlap between F121-9 NADs and ciLAD
ol.ciLAD.non_XL.MEF_LAD <- findOverlapsOfPeaks(bed.in[c(1,3)])
Aizhan.peaks.in.unique <- getPeaksInUnique(ol.ciLAD.non_XL.MEF_LAD)
Aizhan.peaks.in.overlap <- getPeaksInOverlap(ol.ciLAD.non_XL.MEF_LAD,1)
Aizhan.3.sets.of.peaks <- c(Aizhan.peaks.in.unique,Aizhan.peaks.in.overlap[[3]])
names(Aizhan.3.sets.of.peaks)[3] <- names(Aizhan.peaks.in.overlap)[3]
output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4Aizhan_F121_9_avesig_filtered_NADsVenn"
outGrl(Aizhan.3.sets.of.peaks,speci="Mm",output.file.dir)

# 1) Overlap between F121-9 NADs and MEF XL NADs
# 4) MEF XL NADs not shared with F121-9 NADs
# 5) F121-9 NADs not shared with MEF XL NADs
ol.ciLAD.non_XL.MEF_LAD <- findOverlapsOfPeaks(bed.in[c(1,4)])
Aizhan.peaks.in.unique <- getPeaksInUnique(ol.ciLAD.non_XL.MEF_LAD)
Aizhan.peaks.in.overlap <- getPeaksInOverlap(ol.ciLAD.non_XL.MEF_LAD,1)
Aizhan.3.sets.of.peaks <- c(Aizhan.peaks.in.unique,Aizhan.peaks.in.overlap[[3]])
names(Aizhan.3.sets.of.peaks)[3] <- names(Aizhan.peaks.in.overlap)[3]
output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4Aizhan_F121_9_avesig_filtered_NADsVenn"
outGrl(Aizhan.3.sets.of.peaks,speci="Mm",output.file.dir)

# intersectBED(bed.in,1,4)
# intersectBED(bed.in,1,2)
# intersectBED(bed.in,1,3)

intersectBED <- function(bed.in,index1,index2) {
  peaks1 <- bed.in[[index1]]
  peaks2 <- bed.in[[index2]]
  peaks1.and.peaks2 <- GenomicRanges::intersect(peaks1,peaks2)
  outputFileNamePrefix = "common"
  NameX <- paste(names(bed.in)[c(index1,index2)],collapse = "-and-")
  NameY <- paste(c(outputFileNamePrefix,NameX),collapse = "-in-")
  re <- list(peaks1.and.peaks2=peaks1.and.peaks2)
  names(re) <- NameY
  output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4Aizhan_F121_9_avesig_filtered_NADsVenn"
  outGrl(re,speci="Mm",output.file.dir)
}

#' uniqueBED
#'
#' @param bed.in 
#' @param index1
#' @param index2 
#' @param output.file.dir 
#'
#' @return a list including unique peaks
#' 
#' @export 
#'
#' @examples
#' unique.bed <- uniqueBED(bed.in,1,4,output.file.dir)
#' 
uniqueBED <- function(bed.in,index1,index2,output.file.dir) {
  
  peaks1 <- bed.in[[index1]]
  peaks2 <- bed.in[[index2]]
  
  peaks1.only <- GenomicRanges::setdiff(peaks1,peaks2)
  peaks2.only <- GenomicRanges::setdiff(peaks2,peaks1)
  
  grl <- list(peaks1.only=peaks1.only,peaks2.only=peaks2.only)
  names(grl) <- paste0("unique-",names(bed.in)[c(index1,index2)])
  
  outGrl(grl,speci="Mm",output.file.dir)
  grl
}

# # 1) F121-9 vs cLAD vs ciLAD analysis:
#   
# Dropbox (UMass Medical School)\pub_old_galaxy\F121-9\overlap analysis\F121-9 vs CLAD vs ciLAD\common-in-F121_9_and_ciLAD.bed

input.file5 <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/overlap analysis/F121-9 vs CLAD vs ciLAD/common-in-F121_9_and_ciLAD.bed"

# Dropbox (UMass Medical School)\pub_old_galaxy\F121-9\overlap analysis\F121-9 vs CLAD vs ciLAD\common-in-F121_9_and_cLAD.bed

input.file6 <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/overlap analysis/F121-9 vs CLAD vs ciLAD/common-in-F121_9_and_cLAD.bed"

# 2) F121-9 vs. MEF XL NADs analysis:
#   
# Dropbox (UMass Medical School)\pub_old_galaxy\F121-9\overlap analysis\F121-9 vs MEF XL NADs\common-in-F121_9_and-XL_MEF.bed
 
input.file7 <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/overlap analysis/F121-9 vs MEF XL NADs/common-in-F121_9_and-XL_MEF.bed"

# Dropbox (UMass Medical School)\pub_old_galaxy\F121-9\overlap analysis\F121-9 vs MEF XL NADs\unique-F121_9.bed

input.file8 <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/overlap analysis/F121-9 vs MEF XL NADs/unique-F121_9.bed"

# Dropbox (UMass Medical School)\pub_old_galaxy\F121-9\overlap analysis\F121-9 vs MEF XL NADs\unique-XL_MEF.bed

input.file9 <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/overlap analysis/F121-9 vs MEF XL NADs/unique-XL_MEF.bed"

input.files.greater.than.50kb <- c(input.file5,input.file6,input.file7,input.file8,input.file9)

bed.in.50kb.above<-lapply(input.files.greater.than.50kb,function(u){
  peaks=read.table(u,header=F,skip=1)
  colnames(peaks)[1:3]= c("chr","start","end")
  peaks=toGRanges(peaks)
  peaks
})

basename(input.files.greater.than.50kb)

names(bed.in.50kb.above) <- basename(input.files.greater.than.50kb)

names(bed.in.50kb.above)

peak.index <- c(1,2)
name <- tools::file_path_sans_ext(names(bed.in.50kb.above)[c(1,2)])
output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4Aizhan_F121_9_avesig_filtered_NADsVenn"
getCount4Venn(bed.in.50kb.above,peak.index,name,output.file.dir)


peak.index <- c(3,4,5)
name <- tools::file_path_sans_ext(names(bed.in.50kb.above)[c(3,4,5)])
output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4Aizhan_F121_9_avesig_filtered_NADsVenn"
getCount4Venn(bed.in.50kb.above,peak.index,name,output.file.dir)

output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4Aizhan_F121_9_avesig_filtered_NADsVenn/filteredByG50KB"

xx <- lapply(bed.in, function(u){
  x <- u[which(width(u)>50000)]
})

peak.index <- c(1,2,3,4)
name <- c("F121_9","cLAD","ciLAD","XL_MEF")
getCount4Venn(xx,peak.index,name,output.file.dir)

peak.index <- c(1,4)
name <- c("F121_9","XL_MEF")
getCount4Venn(xx,peak.index,name,output.file.dir)

peak.index <- c(1,2,3)
name <- c("F121_9","cLAD","ciLAD")
getCount4Venn(xx,peak.index,name,output.file.dir)

#nonXL.3.subsets <- getUniquePeaks(ol.ciLAD.non_XL.MEF_LAD)

# Hi Aimin,
# 
# Thank you again for the Venn diagrams.
# 
# You have done NAD vs. RNAseq fpkm box plot analysis before (I've attached it to this email). Can you please re-do the box plot, since our F121-9 NAD files have changed after using AveSig>1.7 cutoff. We need 2 box plots, one is NADs vs. MEF RNAseq fpkm, the other is NADs vs. F121-9 RNAseq fpkm values. I incuded all needed file paths below. Please let me know if you have questions.
#                                                             
#                                                             1) For MEF box plot, please use MEF RNAseq from Collas lab (GSM1621026):
#                                                             Dropbox (UMass Medical School)\NADfinder\Aimin\F121_9_RNA_Seq_data\GSM1621026_WT_fpkm

GSM1621026.WT.fpkm <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/F121_9_RNA_Seq_data/GSM1621026_WT_fpkm.csv"

#                                                             
#                                                             F121-9 only NADs:
#                                                             Dropbox (UMass Medical School)\pub_old_galaxy\F121-9\overlap analysis\F121-9 vs MEF XL NADs\unique-F121_9

unique.F121_9 <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/overlap analysis/F121-9 vs MEF XL NADs/unique-F121_9.bed"

Dropbox (UMass Medical School)\pub_old_galaxy\F121-9\overlap analysis\F121-9 vs MEF XL NADs\unique-F121_9

#                                                             
#                                                             MEF XL only NADs:
#                                                             Dropbox (UMass Medical School)\pub_old_galaxy\F121-9\overlap analysis\F121-9 vs MEF XL NADs\unique-XL_MEF

unique.XL_MEF <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/overlap analysis/F121-9 vs MEF XL NADs/unique-XL_MEF.bed"

Dropbox (UMass Medical School)\pub_old_galaxy\F121-9\overlap analysis\F121-9 vs MEF XL NADs\unique-XL_MEF

#                                                             
#                                                             F121-9 and MEF XL shared NADs:
#                                                             Dropbox (UMass Medical School)\pub_old_galaxy\F121-9\overlap analysis\F121-9 vs MEF XL NADs\common-in-F121_9_and-XL_MEF

common.in.F121_9_and_XL_MEF <- "/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/F121-9/overlap analysis/F121-9 vs MEF XL NADs/common-in-F121_9_and-XL_MEF.bed"

Dropbox (UMass Medical School)\pub_old_galaxy\F121-9\overlap analysis\F121-9 vs MEF XL NADs\common-in-F121_9_and-XL_MEF

# NADs (for MEF box plot, please use MEF XL NADs, path shown below):
#   Dropbox (UMass Medical School)\NADfinder\BedFiles\Aizhan F121-9 comparisons\XL_MEF_1.5FC_minp_countF200_auto_q_less_than_e-3_chrX_q_less_than_e-2_adjPval copy _greater_than_50_kb

MEF.XL.NAD <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/BedFiles/Aizhan F121-9 comparisons/XL_MEF_1.5FC_minp_countF200_auto_q_less_than_e-3_chrX_q_less_than_e-2_adjPval copy _greater_than_50_kb.txt"

Dropbox (UMass Medical School)\NADfinder\BedFiles\Aizhan F121-9 comparisons\XL_MEF_1.5FC_minp_countF200_auto_q_less_than_e-3_chrX_q_less_than_e-2_adjPval copy _greater_than_50_kb

#                                                             
#                                                             nonNAD: I don’t have separate files for the “nonNAD” genome subsets you had previously calculated. Can you please calculate nonNADs according to MEF XL NADs:

#                                                             Dropbox (UMass Medical School)\NADfinder\BedFiles\Aizhan F121-9 comparisons\XL_MEF_1.5FC_minp_countF200_auto_q_less_than_e-3_chrX_q_less_than_e-2_adjPval copy _greater_than_50_kb

#                                                             
#                                                             2) For F121-9 box plot, please use F121-9 RNAseq. 
#                                                             Previously, I've sent you bam files for two replicates of F121-9 RNAseq, they're here:
#                                                             Dropbox (UMass Medical School)\pub_old_galaxy\F121-9\Aizhan.allele.specific\RNASEQ\hybrid_STAR\repl_1.bam
#                                                             
#                                                             Dropbox (UMass Medical School)\pub_old_galaxy\F121-9\Aizhan.allele.specific\RNASEQ\hybrid_STAR\repl_2.bam
#                                                             
#                                                             You wrote us that: "I used 2 BAM files from Aizhan to get raw counts, then calculated the FPKM for each replicate, and took the average FPKM to represent gene FPKM value". I don't know where the file for your calculated fpkm values is located. 
#                                                             
#                                                             F121-9 only NADs:
#                                                               Dropbox (UMass Medical School)\pub_old_galaxy\F121-9\overlap analysis\F121-9 vs MEF XL NADs\unique-F121_9

#                                                             
#                                                             MEF XL only NADs:
#                                                               Dropbox (UMass Medical School)\pub_old_galaxy\F121-9\overlap analysis\F121-9 vs MEF XL NADs\unique-XL_MEF
#                                                             
#                                                             F121-9 and MEF XL shared NADs:
#                                                               Dropbox (UMass Medical School)\pub_old_galaxy\F121-9\overlap analysis\F121-9 vs MEF XL NADs\common-in-F121_9_and-XL_MEF
#                                                             
#                                                             NADs (for F121-9 box plot, please use F121-9 NADs, path shown below):
#                                                               Dropbox (UMass Medical School)\pub_old_galaxy\F121-9\avesig filtered NADs\F121_9_all_1.5FC_F200_minp_adjpVal_avesig_greater_than_1.7_NAD_greater_50kb

F121_9.NAD <- "~/Aimin/DropboxUmass/pub_old_galaxy/F121-9/avesig filtered NADs/F121_9_all_1.5FC_F200_minp_adjpVal_avesig_greater_than_1.7_NAD_greater_50kb.bed"

#                                                             
#                                                             nonNAD: please calculate nonNADs according to F121-9 NADs:
#                                                               Dropbox (UMass Medical School)\pub_old_galaxy\F121-9\avesig filtered NADs\F121_9_all_1.5FC_F200_minp_adjpVal_avesig_greater_than_1.7_NAD_greater_50kb
#                                                             
#                                                             Thank you!
#                                                               Aizhan
#                                                             

input.files.new <- c(unique.F121_9,unique.XL_MEF,common.in.F121_9_and_XL_MEF,MEF.XL.NAD,F121_9.NAD)
bed.in.new<-lapply(1:length(input.files.new),function(u,input.files.new){

  x <- input.files.new[[u]]
  
  if(u < 4){peaks=read.table(x,header=F)}else{peaks=read.table(x,header=F,skip=1)}
  
  colnames(peaks)[1:3]= c("chr","start","end")
  peaks=toGRanges(peaks)
  peaks
},input.files.new)

names(bed.in.new) <- basename(input.files.new)
names(bed.in.new) <- c("unique_F121_9","unique_XL_MEF","F121_9_and_XL_MEF","XL_NAD","F121_9_NAD")

output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4Aizhan_F121_9_avesig_filtered_NADsVenn/MEF_box_plot"

bed.in.new.anno <- getAnnotatedGene(bed.in.new[c(1:3,4)],"Mm")
bed.in.new.anno.YYY <- getFPKM4DiffSet(bed.in.new.anno,fpkm.value)
getBoxPlot4FPKMOfSubsetPeaks7(bed.in.new.anno.YYY,fpkm.value,"XL_NAD","GSM1621026_WT_FPKM.txt",output.file.dir)
save.image(file = file.path(output.file.dir,paste0(basename(output.file.dir),".RData")))
savehistory(file = file.path(output.file.dir,paste0(basename(output.file.dir),".Rhistory")))


output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4Aizhan_F121_9_avesig_filtered_NADsVenn/F121-9_box_plot"
bed.in.new.anno <- getAnnotatedGene(bed.in.new[c(1:3,5)],"Mm")
bed.in.new.anno.YYY <- getFPKM4DiffSet(bed.in.new.anno,rna.seq.data.from.Aizhan)
getBoxPlot4FPKMOfSubsetPeaks7(bed.in.new.anno.YYY,rna.seq.data.from.Aizhan,"F121_9_NAD","rnaSeq_data_for_F121_9_FPKM.txt",output.file.dir)
save.image(file = file.path(output.file.dir,paste0(basename(output.file.dir),".RData")))
savehistory(file = file.path(output.file.dir,paste0(basename(output.file.dir),".Rhistory")))



getBoxPlot4FPKMOfSubsetPeaks2 <- function(YYY,fpkm.value,output.file.dir){
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  fpkm.value.ref <- data.frame(SetName=rep("wholeGenome",dim(fpkm.value)[1]),GeneName=fpkm.value$gene.name,FPKM=fpkm.value$fpkm)
  
  NAD <- bed.in.new.anno.YYY[which(bed.in.new.anno.YYY$SetName=="NAD"),]
  
  NAD.fpkm.value <- fpkm.value.ref[which(fpkm.value.ref$GeneName %in% NAD$GeneName),]
  NAD.fpkm.value$SetName <- "NAD"
  nonNAD.fpkm.value <- fpkm.value.ref[-which(fpkm.value.ref$GeneName %in% YYY$GeneName),]
  nonNAD.fpkm.value$SetName <- "nonNAD"
  
  YYYY <- rbind(YYY,NAD.fpkm.value,nonNAD.fpkm.value,fpkm.value.ref)
  
  png(file = file.path(output.file.dir,paste0(paste(levels(YYYY$SetName),collapse = "-"),".png")),width = 1500, height = 480)
  boxplot(log10(as.numeric(as.character(YYYY$FPKM))+1)~YYYY$SetName,ylab = "log10(FPKM+1)")
  dev.off()
  
  YYYY
  
}

output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4Aizhan_F121_9_avesig_filtered_NADsVenn/MEF_box_plot_linear"
bed.in.new.anno.YYY <- getFPKM4DiffSet(bed.in.new.anno,fpkm.value)
getBoxPlot4FPKMOfSubsetPeaks3(bed.in.new.anno.YYY,fpkm.value,output.file.dir)

output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4Aizhan_F121_9_avesig_filtered_NADsVenn/F121-9_box_plot_linear"
bed.in.new.anno.YYY <- getFPKM4DiffSet(bed.in.new.anno,rna.seq.data.from.Aizhan)
getBoxPlot4FPKMOfSubsetPeaks3(bed.in.new.anno.YYY,rna.seq.data.from.Aizhan,output.file.dir)

getBoxPlot4FPKMOfSubsetPeaks3 <- function(YYY,fpkm.value,output.file.dir){
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  fpkm.value.ref <- data.frame(SetName=rep("wholeGenome",dim(fpkm.value)[1]),GeneName=fpkm.value$gene.name,FPKM=fpkm.value$fpkm)
  
  NAD <- bed.in.new.anno.YYY[which(bed.in.new.anno.YYY$SetName=="NAD"),]
  
  NAD.fpkm.value <- fpkm.value.ref[which(fpkm.value.ref$GeneName %in% NAD$GeneName),]
  NAD.fpkm.value$SetName <- "NAD"
  nonNAD.fpkm.value <- fpkm.value.ref[-which(fpkm.value.ref$GeneName %in% YYY$GeneName),]
  nonNAD.fpkm.value$SetName <- "nonNAD"
  
  YYYY <- rbind(YYY,NAD.fpkm.value,nonNAD.fpkm.value,fpkm.value.ref)
  
  png(file = file.path(output.file.dir,paste0(paste(levels(YYYY$SetName),collapse = "-"),".png")),width = 1500, height = 480)
  boxplot(as.numeric(as.character(YYYY$FPKM))~YYYY$SetName,ylab = "FPKM")
  dev.off()

YYYY

}

output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4Aizhan_F121_9_avesig_filtered_NADsVenn/MEF_box_plot_linear_rm_gt_100"
bed.in.new.anno.YYY <- getFPKM4DiffSet(bed.in.new.anno,fpkm.value)
getBoxPlot4FPKMOfSubsetPeaks4(bed.in.new.anno.YYY,fpkm.value,output.file.dir)

output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4Aizhan_F121_9_avesig_filtered_NADsVenn/F121-9_box_plot_linear_rm_gt_100"
bed.in.new.anno.YYY <- getFPKM4DiffSet(bed.in.new.anno,rna.seq.data.from.Aizhan)
getBoxPlot4FPKMOfSubsetPeaks4(bed.in.new.anno.YYY,rna.seq.data.from.Aizhan,output.file.dir)

getBoxPlot4FPKMOfSubsetPeaks4 <- function(YYY,fpkm.value,output.file.dir){
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  fpkm.value.ref <- data.frame(SetName=rep("wholeGenome",dim(fpkm.value)[1]),GeneName=fpkm.value$gene.name,FPKM=fpkm.value$fpkm)
  
  NAD <- bed.in.new.anno.YYY[which(bed.in.new.anno.YYY$SetName=="NAD"),]
  
  NAD.fpkm.value <- fpkm.value.ref[which(fpkm.value.ref$GeneName %in% NAD$GeneName),]
  NAD.fpkm.value$SetName <- "NAD"
  nonNAD.fpkm.value <- fpkm.value.ref[-which(fpkm.value.ref$GeneName %in% YYY$GeneName),]
  nonNAD.fpkm.value$SetName <- "nonNAD"
  
  YYYY <- rbind(YYY,NAD.fpkm.value,nonNAD.fpkm.value,fpkm.value.ref)
  
  png(file = file.path(output.file.dir,paste0(paste(levels(YYYY$SetName),collapse = "-"),".png")),width = 1500, height = 480)
  
  YYYY <- YYYY[which(as.numeric(as.character(YYYY$FPKM))<=100),]
  
  boxplot(as.numeric(as.character(YYYY$FPKM))~YYYY$SetName,ylab = "FPKM")
  
  dev.off()
  
  YYYY
  
}

output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4Aizhan_F121_9_avesig_filtered_NADsVenn/MEF_box_plot_linear_rm_gt_100_log"
bed.in.new.anno.YYY <- getFPKM4DiffSet(bed.in.new.anno,fpkm.value)
getBoxPlot4FPKMOfSubsetPeaks5(bed.in.new.anno.YYY,fpkm.value,output.file.dir)

output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4Aizhan_F121_9_avesig_filtered_NADsVenn/F121-9_box_plot_linear_rm_gt_100_log"
bed.in.new.anno.YYY <- getFPKM4DiffSet(bed.in.new.anno,rna.seq.data.from.Aizhan)
getBoxPlot4FPKMOfSubsetPeaks5(bed.in.new.anno.YYY,rna.seq.data.from.Aizhan,output.file.dir)

getBoxPlot4FPKMOfSubsetPeaks5 <- function(YYY,fpkm.value,output.file.dir){
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  fpkm.value.ref <- data.frame(SetName=rep("wholeGenome",dim(fpkm.value)[1]),GeneName=fpkm.value$gene.name,FPKM=fpkm.value$fpkm)
  
  NAD <- bed.in.new.anno.YYY[which(bed.in.new.anno.YYY$SetName=="NAD"),]
  
  NAD.fpkm.value <- fpkm.value.ref[which(fpkm.value.ref$GeneName %in% NAD$GeneName),]
  NAD.fpkm.value$SetName <- "NAD"
  nonNAD.fpkm.value <- fpkm.value.ref[-which(fpkm.value.ref$GeneName %in% YYY$GeneName),]
  nonNAD.fpkm.value$SetName <- "nonNAD"
  
  YYYY <- rbind(YYY,NAD.fpkm.value,nonNAD.fpkm.value,fpkm.value.ref)
  
  png(file = file.path(output.file.dir,paste0(paste(levels(YYYY$SetName),collapse = "-"),".png")),width = 1500, height = 480)
  
  YYYY <- YYYY[which(as.numeric(as.character(YYYY$FPKM))<=100),]
  
  boxplot(log10(as.numeric(as.character(YYYY$FPKM))+1)~YYYY$SetName,ylab = "log10(FPKM+1)")
  
  dev.off()
  
  YYYY
  
}

GSE90894.RPKM <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/F121_9_RNA_Seq_data/GSE90894_RPKM_mRNAseq_table.xlsx"
library(readxl)
GSE90894.RPKM.data <- read_xlsx(GSE90894.RPKM,sheet=1)
colnames(GSE90894.RPKM.data) <- GSE90894.RPKM.data[4,]

GSE90894.RPKM.data.1 <- GSE90894.RPKM.data[-c(1:4),]
colnames(GSE90894.RPKM.data.1)

rpkm.MEFs <- as.data.frame(GSE90894.RPKM.data.1[,c(1,2)])
rpkm.ESCs <- as.data.frame(GSE90894.RPKM.data.1[,c(1,6)])
colnames(rpkm.MEFs) <- c("gene.name","fpkm")
colnames(rpkm.ESCs) <- c("gene.name","fpkm")

countToFpkm <- function(counts, effLen)
{
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

input.count.file <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/F121_9_RNA_Seq_data/GSM1418813_ES1_counts.txt"
  
GSM1418813.ES1.fpkm <- countFile2Fpkm(input.count.file)

countFile2Fpkm <- function(input.count.file) {
  
  GSM1418813.ES1.counts <- input.count.file
  count.data <- read.table(GSM1418813.ES1.counts,header= F)
  colnames(count.data) <- c("gene.name","count")
  mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
  G_list <- getBM(filters= "mgi_symbol", attributes= c("start_position","end_position","mgi_symbol"),values=count.data$gene.name,mart= mart)
  geneNameLen <- data.frame(gene.name=G_list$mgi_symbol,len=G_list$end_position-G_list$start_position)
  geneNameCountLen <- merge(count.data,geneNameLen,by="gene.name")
  GSM1418813.ES1.fpkm <- data.frame(gene.name=geneNameCountLen$gene.name,fpkm=countToFpkm(geneNameCountLen[,2],geneNameCountLen[,3]))
  GSM1418813.ES1.fpkm
  
}

output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4Aizhan_F121_9_avesig_filtered_NADsVenn/MEF_Plath_GSE90984_box_plot"

bed.in.new.anno <- getAnnotatedGene(bed.in.new[c(1:3,4)],"Mm")
bed.in.new.anno.YYY <- getFPKM4DiffSet(bed.in.new.anno,rpkm.MEFs)
getBoxPlot4FPKMOfSubsetPeaks7(bed.in.new.anno.YYY,rpkm.MEFs,"XL_NAD","rpkm.MEFs.txt",output.file.dir)
save.image(file = file.path(output.file.dir,paste0(basename(output.file.dir),".RData")))
savehistory(file = file.path(output.file.dir,paste0(basename(output.file.dir),".Rhistory")))

#bed.in.new.anno.YYY <- getFPKM4DiffSet(bed.in.new.anno,rpkm.MEFs)
#getBoxPlot4FPKMOfSubsetPeaks6(bed.in.new.anno.YYY,rpkm.MEFs,"XL_NAD",output.file.dir)


output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4Aizhan_F121_9_avesig_filtered_NADsVenn/mESC_Plath_GSE90984_box_plot"

bed.in.new.anno <- getAnnotatedGene(bed.in.new[c(1:3,5)],"Mm")
bed.in.new.anno.YYY <- getFPKM4DiffSet(bed.in.new.anno,rpkm.ESCs)
getBoxPlot4FPKMOfSubsetPeaks7(bed.in.new.anno.YYY,rpkm.ESCs,"F121_9_NAD","rpkm.ESCs.txt",output.file.dir)
save.image(file = file.path(output.file.dir,paste0(basename(output.file.dir),".RData")))
savehistory(file = file.path(output.file.dir,paste0(basename(output.file.dir),".Rhistory")))

#bed.in.new.anno.YYY <- getFPKM4DiffSet(bed.in.new.anno,rpkm.ESCs)
#getBoxPlot4FPKMOfSubsetPeaks6(bed.in.new.anno.YYY,rpkm.ESCs,"F121_9_NAD",output.file.dir)

output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4Aizhan_F121_9_avesig_filtered_NADsVenn/mESC_Lowe_GSM1418813_box_plot"

bed.in.new.anno <- getAnnotatedGene(bed.in.new[c(1:3,5)],"Mm")
bed.in.new.anno.YYY <- getFPKM4DiffSet(bed.in.new.anno,GSM1418813.ES1.fpkm)
getBoxPlot4FPKMOfSubsetPeaks7(bed.in.new.anno.YYY,GSM1418813.ES1.fpkm,"F121_9_NAD","GSM1418813.ES1.fpkm.txt",output.file.dir)
save.image(file = file.path(output.file.dir,paste0(basename(output.file.dir),".RData")))
savehistory(file = file.path(output.file.dir,paste0(basename(output.file.dir),".Rhistory")))

#bed.in.new.anno.YYY <- getFPKM4DiffSet(bed.in.new.anno,GSM1418813.ES1.fpkm)
#getBoxPlot4FPKMOfSubsetPeaks6(bed.in.new.anno.YYY,GSM1418813.ES1.fpkm,"F121_9_NAD",output.file.dir)



getBoxPlot4FPKMOfSubsetPeaks6 <- function(YYY,fpkm.value,setName,output.file.dir){
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  fpkm.value.ref <- data.frame(SetName=rep("wholeGenome",dim(fpkm.value)[1]),GeneName=fpkm.value$gene.name,FPKM=fpkm.value$fpkm)
  
  NAD <- YYY[which(YYY$SetName==setName),]
  
  NAD.fpkm.value <- fpkm.value.ref[which(fpkm.value.ref$GeneName %in% NAD$GeneName),]
  NAD.fpkm.value$SetName <- "NAD"
  nonNAD.fpkm.value <- fpkm.value.ref[-which(fpkm.value.ref$GeneName %in% YYY$GeneName),]
  nonNAD.fpkm.value$SetName <- "nonNAD"
  
  YYYY <- rbind(YYY,NAD.fpkm.value,nonNAD.fpkm.value,fpkm.value.ref)
  
  png(file = file.path(output.file.dir,paste0(paste(levels(YYYY$SetName),collapse = "-"),".png")),width = 1500, height = 480)
  
  #YYYY <- YYYY[which(as.numeric(as.character(YYYY$FPKM))<=100),]
  
  boxplot(log10(as.numeric(as.character(YYYY$FPKM))+1)~YYYY$SetName,ylab = "log10(FPKM+1)")
  
  dev.off()
  
  YYYY
  
}

input.files.6.4.2019 <- c(input.file6,input.file5,input.file1)
gRName <-  c("unique_F121_9","unique_XL_MEF","F121_9_and_XL_MEF","XL_NAD","F121_9_NAD")

bedFiless2Gr <- function(input.files.new,gRName) {
  
  bed.in.new<-lapply(1:length(input.files.new),function(u,input.files.new){
    
    x <- input.files.new[[u]]
    
    if(u <=2 ){peaks=read.table(x,header=F)}else{peaks=read.table(x,header=F,skip=1)}
    
    colnames(peaks)[1:3]= c("chr","start","end")
    peaks=toGRanges(peaks)
    peaks
  },input.files.new)
  
  names(bed.in.new) <- basename(input.files.new)
  names(bed.in.new) <- gRName
  bed.in.new
}

gRName <- c("F121_9_Type_I_NADs","F121_9_Type_II_NADs","F121_9_NADs")
gr.6.4.2019 <- bedFiless2Gr(input.files.6.4.2019,gRName) 

gr.6.4.2019.anno <- getAnnotatedGene(gr.6.4.2019,"Mm")

output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output4_F121_9_Type_I_NADs_F121_9_Type_II_NADs_F121_9_NADs_F121_9_nonNADs_box_plot"

gr.6.4.2019.anno.YYY <- getFPKM4DiffSet(gr.6.4.2019.anno,rna.seq.data.from.Aizhan)
getBoxPlot4FPKMOfSubsetPeaks7(gr.6.4.2019.anno.YYY,rna.seq.data.from.Aizhan,"F121_9_NADs","rnaSeq_data_from_Aizhan_FPKM.txt",output.file.dir)

getBoxPlot4FPKMOfSubsetPeaks7 <- function(YYY,fpkm.value,setName,rnaSeqFPKMFileName,output.file.dir){
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  write.table(fpkm.value,file=file.path(output.file.dir,rnaSeqFPKMFileName),sep="\t",quote = FALSE,row.names = FALSE,col.names = TRUE)
  
  fpkm.value.ref <- data.frame(SetName=rep("wholeGenome",dim(fpkm.value)[1]),GeneName=fpkm.value$gene.name,FPKM=fpkm.value$fpkm)
  
  #setName <- "F121_9_NADs"
  NAD <- YYY[which(YYY$SetName==setName),]
  
  #NAD.fpkm.value <- fpkm.value.ref[which(fpkm.value.ref$GeneName %in% NAD$GeneName),]
  #NAD.fpkm.value$SetName <- "NAD"
  nonNAD.fpkm.value <- fpkm.value.ref[-which(fpkm.value.ref$GeneName %in% YYY$GeneName),]
  nonNAD.fpkm.value$SetName <- "nonNAD"
  
  YYYY <- rbind(YYY,NAD,nonNAD.fpkm.value,fpkm.value.ref)
  
  png(file = file.path(output.file.dir,paste0(paste(levels(YYYY$SetName),collapse = "-"),".png")),width = 1500, height = 480)
  
  boxplot(log10(as.numeric(as.character(YYYY$FPKM))+1)~YYYY$SetName,ylab = "log10(FPKM+1)")
  
  dev.off()
  
  write.table(YYYY,file=file.path(output.file.dir,"SetName_FPKM.txt"),sep="\t",quote = FALSE,row.names = FALSE,col.names = TRUE)
  
  YYYY
  
}

gr.whole.genome <- GRanges(seqinfo(BSgenome.Mmusculus.UCSC.mm10))
F121_9_nonNADs <- list(F121_9_nonNADs=GenomicRanges::setdiff(gr.whole.genome,gr.6.4.2019[[3]]))
outGrl(F121_9_nonNADs,speci="Mm",output.file.dir)

save.image(file = file.path(output.file.dir,"6_4_2019_FPKM.RData"))
savehistory(file = file.path(output.file.dir,"6_4_2019_FPKM.Rhistory"))

