#!/usr/bin/env Rscript

#Here's a file with F121-9 NADs and AveSig values for each peak:
#Dropbox (UMass Medical School)\NADfinder\tileCount-NADfinder1.5.2\backgroundCorrected\F121_9adjPlt0.05-1.5FC-countF200-minimump.xls

F121_9.NAD.AveSig <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/tileCount-NADfinder1.5.2/backgroundCorrected/F121_9adjPlt0.05-1.5FC-countF200-minimump-2.xls"

#And here's the file with crosslinked MEF NADs and AveSig values for each peak:

MEF.NAS.AveSig <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/tileCount-NADfinder1.5.2/backgroundCorrected/CrossLinkedMEFadjPlt0.05-1.5FC-minimump-countF200-2.xls"

data.F121_9.NAD.AveSig <- read.xls(F121_9.NAD.AveSig)
data.MEF.NAS.AveSig <- read.xls(MEF.NAS.AveSig)

F121_9.and.MEF <- GenomicRanges::intersect(toGRanges(data.F121_9.NAD.AveSig[,c(1:3,8)]),toGRanges(data.MEF.NAS.AveSig[,c(1:3,8)]))

F121_9.and.MEF

hits <- findOverlaps(toGRanges(data.F121_9.NAD.AveSig[,c(1:3)]),F121_9.and.MEF)
data.F121_9.NAD.AveSig.part.in.overlap <- toGRanges(data.F121_9.NAD.AveSig[,c(1:3,8)])[queryHits(hits)]                    
hits <- findOverlaps(toGRanges(data.MEF.NAS.AveSig[,c(1:3)]),F121_9.and.MEF)
data.MEF.NAS.AveSig.part.in.overlap <- toGRanges(data.MEF.NAS.AveSig[,c(1:3,9)])[queryHits(hits)]                 

x <- mcols(data.F121_9.NAD.AveSig.part.in.overlap)
y <- mcols(data.MEF.NAS.AveSig.part.in.overlap)
   
Gr2Df <- function(data.F121_9.NAD.AveSig.part.in.overlap) {
  data.F121_9.NAD.AveSig.part.in.overlap.1 <- data.frame(seqnames=seqnames(data.F121_9.NAD.AveSig.part.in.overlap),
                                                         start=start(data.F121_9.NAD.AveSig.part.in.overlap),end=end(data.F121_9.NAD.AveSig.part.in.overlap),AveSig=mcols(data.F121_9.NAD.AveSig.part.in.overlap)[,1])
}

data.F121_9.NAD.AveSig.part.in.overlap.1 <- Gr2Df(data.F121_9.NAD.AveSig.part.in.overlap)
data.MEF.NAS.AveSig.part.in.overlap.1 <- Gr2Df(data.MEF.NAS.AveSig.part.in.overlap)

write.table(data.F121_9.NAD.AveSig.part.in.overlap.1,file = file.path(output.dir,paste0("F121_9_NAD_AveSig_part_in_overlap",".bed")),append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = F)

write.table(data.MEF.NAS.AveSig.part.in.overlap.1,file = file.path(output.dir,paste0("MEF_NAD_AveSig_part_in_overlap",".bed")),append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = F)

output.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/tileCount-NADfinder1.5.2/backgroundCorrected/Aimin/OutPut"
if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}

pdf(file.path(output.dir,"Comparing_F121_9_NADs_with_MEF_NADs.pdf"))
plot(y[,1], x[,1], main="Comparing F121-9 NADs with MEF NADs", 
     xlab="MEF_NAD_AveSig", ylab="F121-9_NAD_AveSig", pch=19)
dev.off()

save.image(file = file.path(output.dir,"Comparing_F121_9_NADs_with_MEF_NADs.RData"))
savehistory(file = file.path(output.dir,"Comparing_F121_9_NADs_with_MEF_NADs.Rhistory"))
