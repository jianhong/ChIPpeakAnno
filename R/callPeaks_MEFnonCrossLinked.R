library(NADfinder)
#setwd("/project/umw_paul_kaufman/jz57w/NADfinder_1.5.1_peaks")
setwd("~/Dropbox (UMass Medical School)/NADfinder/tileCount-NADfinder1.5.2")
#data(triplicate.count)
#se <- triplicate.count

#Non-crosslinked genomic inputs (samples #18, 29, 33):


se18 <- readRDS("1.tileCounts.S_18.RDS")
se29 <- readRDS("1.tileCounts.S_29.RDS")
se33 <- readRDS("1.tileCounts.S_33.RDS")
 
se <- cbind(se18, se29, se33)
metadata(se)$lib.size.chrom <- cbind(metadata(se18)$lib.size.chrom, 
				  metadata(se29)$lib.size.chrom,
				 metadata(se33)$lib.size.chrom)

nuc.cols <- colnames(se)[c(2,4,6)]
gen.cols <- colnames(se)[c(1,3,5)]
## Calculate ratios for nucleoloar vs genomic samples.

se.log2 <- log2se(se, 
             nucleolusCols = nuc.cols,
             genomeCols = gen.cols,
             transformation = "log2CPMRatio")

se.m <- as.data.frame(assays(se))

se.log2.0 <- log2se(se,nucleolusCols = nuc.cols,genomeCols = gen.cols,transformation = "log2Ratio")

se.log2.m <- as.data.frame(assays(log2se(se, 
                                         nucleolusCols = nuc.cols,
                                         genomeCols = gen.cols,
                                         transformation = "log2Ratio")))


se.log2.1 <- as.data.frame(assays(log2se(se, 
                  nucleolusCols = nuc.cols,
                  genomeCols = gen.cols,
                  transformation = "log2Ratio")))

seList<- smoothRatiosByChromosome(se, chrom.level.background = FALSE)
saveRDS(seList,file= file.path("~/Aimin/DropboxUmass/NADfinder/Aimin/BoundaryAnalysis","seLits.rds"))

test = 0

if (test)
{
se18 <- log2se(se18, nucleolusCols = colnames(se18)[2],
             genomeCols = colnames(se18)[1],
             transformation = "log2CPMRatio")

se29 <- log2se(se29, nucleolusCols = colnames(se29)[2],
             genomeCols = colnames(se29)[1],
             transformation = "log2CPMRatio")

se33 <- log2se(se33, nucleolusCols = colnames(se33)[2],
             genomeCols = colnames(se33)[1],
             transformation = "log2CPMRatio") 

#### double check the ratio calculation 
if(abs(sum((assays(se29)$ratio - assays(se)$ratio[,2]))) > 0.000000000001) {1}
if(abs(sum((assays(se33)$ratio - assays(se)$ratio[,3]))) > 0.000000000001) {1}
if(abs(sum((assays(se18)$ratio - assays(se)$ratio[,1]))) > 0.000000000001) {1}
}

#seList<- smoothRatiosByChromosome(se)

#source("../../Bioconductor/Trunk/NADfinder/R/smoothRatiosByChromosome2.R")

#seList<- smoothRatiosByChromosome(se, chrom.level.background = TRUE)


test = 0
#cumulativePercentage(seList[["chr18"]])
if (test)
{
    
peaks <- lapply(seList, callPeaks,
    cutoffAdjPvalue=0.05, countFilter=200, combineP.method ="Fishers")


    #seList<- smoothRatiosByChromosome2(se, chrom.level.background = FALSE, backgroundCorrection =  FALSE)

    peaks <- lapply(seList, callPeaks,
        cutoffAdjPvalue=0.05, countFilter=200, combineP.method = "Browns")
#peaks<- unlist(GRangesList(peaks))
        peaks <- do.call(rbind, lapply(peaks, function(.ele) {if(length(.ele) >0) as.data.frame(.ele)}))

                
                write.table(as.data.frame(peaks), file ="nonCrossLinkedMEFadjPlt0.05-1.5FC-browns-countF200.xls", sep ="\t", row.names = FALSE)
                
peaks <- lapply(seList, callPeaks,
                cutoffAdjPvalue=0.05, countFilter=200, combineP.method = "Browns", lfc=log2(1.2))
            peaks <- do.call(rbind, lapply(peaks, function(.ele) {if(length(.ele) >0) as.data.frame(.ele)}))
                
            write.table(as.data.frame(peaks), file ="nonCrossLinkedMEFadjPlt0.05-1.2FC-browns-countF200.xls", sep ="\t", row.names = FALSE)

peaks <- lapply(seList, callPeaks, 
                cutoffAdjPvalue=0.05, countFilter=200, combineP.method = "minimump", lfc=log2(1.5))
        peaks <- do.call(rbind, lapply(peaks, function(.ele) {if(length(.ele) >0) as.data.frame(.ele)}))
    write.table(as.data.frame(peaks), file ="nonCrossLinkedMEFadjPlt0.05-1.5FC-minimump-countF200-noBaseLineCor.xls", sep ="\t", row.names = FALSE)
            
peaks <- lapply(seList, callPeaks,
            cutoffAdjPvalue=0.05, countFilter=200, combineP.method = "minimump", lfc=log2(1.5))
            peaks <- do.call(rbind, lapply(peaks, function(.ele) {if(length(.ele) >0) as.data.frame(.ele)}))
            write.table(as.data.frame(peaks), file ="nonCrossLinkedMEFadjPlt0.05-1.5FC-minimump-countF200.xls", sep ="\t", row.names = FALSE)



########

peaks <- lapply(seList, callPeaks,
                cutoffAdjPvalue=0.05, countFilter=1, combineP.method = "Fishers")

peaks<- unlist(GRangesList(peaks))
write.table(as.data.frame(peaks), file ="nonCrossLinkedMEFadjPlt0.05-1.5FC-Fishers.xls", sep ="\t", row.names = FALSE)

peaks <- lapply(seList, callPeaks,
cutoffAdjPvalue=0.05, countFilter=1, combineP.method = "Fishers", lfc = log2(1.2))

peaks<- unlist(GRangesList(peaks))
write.table(as.data.frame(peaks), file ="nonCrossLinkedMEFadjPlt0.05-1.2FC-Fishers.xls", sep ="\t", row.names = FALSE)


}
peaks <- lapply(seList, callPeaks, 
                cutoffAdjPvalue=0.05, countFilter=1)

peaks<- unlist(GRangesList(peaks))

#write.table(as.data.frame(peaks), file ="nonCrossLinkedMEFadjPlt0.05-minimump.xls", sep ="\t", row.names = FALSE)

write.table(as.data.frame(peaks), file ="nonCrossLinkedMEFadjPlt0.05-meanSigs-genomeLevelBG.xls", sep ="\t", row.names = FALSE)

#peaks.meanSig <- peaks
