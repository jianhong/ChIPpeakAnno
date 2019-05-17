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
