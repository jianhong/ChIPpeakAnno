# OK, let me check which one I used, and I can add  type I or type II NADs to that Figure.
# Best,
# Aimin
# 
# From: "Kaufman, Paul (MolMed-PGF&E)" <Paul.Kaufman1@umassmed.edu>
#   Date: Monday, March 11, 2019 at 1:25 PM
# To: "Yan, Aimin" <Aimin.Yan@umassmed.edu>, "Zhu, Lihua (Julie)" <Julie.Zhu@umassmed.edu>
#   Subject: addition to RT Jaccard analysis
# 
# Hi Aimin- I have a question about the previous jaccard analysis of replication timing (attached). I really like this figure and I think it does a better job than the Venn diagrams we had before. However, I realize now that the figure this was meant to replace also included separate analysis of Type I and Type II NADs. Could we simply add those two lines?  The links to the Type I and II NADs are here (should we use XL or nonXL?- what did you use for NADs for the first batch? I think we can pick just one for now).
# 
# https://www.dropbox.com/s/l702oie52yzyvbz/ciLAD_overlap_XL_MEF_1.5FC_qe-3_chrX_qe-2%20min%20max.bed?dl=0
# https://www.dropbox.com/s/iaxgcr9f29ustyt/ciLAD.overlap-nonXL_MEF_1.5FC_qe-3_chrX_qe-2%20igv.bed?dl=0
# https://www.dropbox.com/s/m2d4gwqij411nfc/XL_MEF_1.5FC_qe-3_chrX_qe-2_overlap-MEFLAD%20min%20max.bed?dl=0
# https://www.dropbox.com/s/xckfkjm2rfns9x1/nonXL_MEF_qe-3_chrX_qe-2.overlap.MEF-LADs%20igv.bed?dl=0
# 
# 
# thanks, paul

# ~/Aimin/DropboxUmass/NADfinder/Aimin/Output/BedFromPaul

generateHMbasedOnDistanceMatrix <- function(re.out.rt.mef.5.sets, output.file.dir, d,m1,m2) {
  
  file.name <- paste(names(re.out.rt.mef.5.sets),collapse = "_")
  
  tiff(file.path(output.file.dir,paste0(file.name,"_heatmap4JC.tiff")))
  colors = c(seq(0,0.001,length=8),seq(0.002,0.008,length=8),seq(0.009,0.2,length=2),seq(0.3,0.4,length=2),seq(0.5,1,length=5))
  
  my_palette <- redgreen(24)
  
  heatmap.3(d,col=my_palette,
            cellnote=round(d,2),
            breaks=colors, density.info="none", trace="none", key=T,symm=T,symkey=F,symbreaks=T,margin=c(m1, m2),key.xtickfun=function() {
              breaks <- parent.frame()$breaks
              return(list(
                at=parent.frame()$scale01(c(breaks[1],breaks[length(breaks)])),
                labels=c(as.character(breaks[1]),as.character(breaks[length(breaks)]))
              ))
            })
  dev.off()
}

re.out.rt.mef.5.sets <- re.out.rt.mef[c(17,16,13,5,14)]
names(re.out.rt.mef.5.sets) <- c("lateRT","earlyRT","nonXL-NAD","ciLAD","LAD")
d <- getJC4aSetOfPeaks(re.out.rt.mef.5.sets, output.file.dir)
generateHMbasedOnDistanceMatrix(re.out.rt.mef.5.sets, output.file.dir,d)

re.out.rt.mef.5.sets <- re.out.rt.mef[c(17,16,15,5,14)]
names(re.out.rt.mef.5.sets) <- c("lateRT","earlyRT","XL-NAD","ciLAD","LAD")
d <- getJC4aSetOfPeaks(re.out.rt.mef.5.sets, output.file.dir)
generateHMbasedOnDistanceMatrix(re.out.rt.mef.5.sets, output.file.dir,d)

re.out.rt.mef.5.sets <- re.out.rt.mef[c(17,16,13,15,5,14)]
names(re.out.rt.mef.5.sets) <- c("lateRT","earlyRT","nonXL-NAD","XL-NAD","ciLAD","LAD")
d <- getJC4aSetOfPeaks(re.out.rt.mef.5.sets, output.file.dir)
generateHMbasedOnDistanceMatrix(re.out.rt.mef.5.sets, output.file.dir,d)

input.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output/BedFromPaul"
bed.in.new.4.JC <- getBedFiles(input.file.dir)
bed.in.new.4.JC.1 <- bed.in.new.4.JC

names(bed.in.new.4.JC) <- c("TypeIIXL","TypeIInonXL","TypeInonXL","TypeIXL")
d <- getJC4aSetOfPeaks(bed.in.new.4.JC, output.file.dir)
generateHMbasedOnDistanceMatrix(bed.in.new.4.JC, output.file.dir,d,10,10)

re.out.rt.mef.10.sets <- c(re.out.rt.mef.5.sets,bed.in.new.4.JC)
d <- getJC4aSetOfPeaks(re.out.rt.mef.10.sets, output.file.dir)
generateHMbasedOnDistanceMatrix(re.out.rt.mef.10.sets, output.file.dir,d,7,7)

d.8.8 <- d[3:10,3:10]
d.2.8 <-d[1:2,3:10] 
d.8.8.u <- unique(c(d.8.8))

d.8.8.u[order(d.8.8.u)]

tiff(file.path(output.file.dir,paste0(file.name,"_heatmap4JC.tiff")))

colors = c(seq(0,0.001,length=8),seq(0.002,0.008,length=8),seq(0.009,0.2,length=2),seq(0.3,0.4,length=2),seq(0.5,1,length=5))

colors <- c(0,d.8.8.u[order(d.8.8.u)])
my_palette <- redgreen(29)
heatmap.3(d.8.8,col=my_palette,
          breaks=colors, density.info="none", trace="none", key=T,symm=T,symkey=F,symbreaks=T,margin=c(9, 9),key.xtickfun=function() {
            breaks <- parent.frame()$breaks
            return(list(
              at=parent.frame()$scale01(c(breaks[1],breaks[length(breaks)])),
              labels=c(as.character(breaks[1]),as.character(breaks[length(breaks)]))
            ))
          })
dev.off()

d.2.8 <-d[1:2,3:10] 
d.2.8.u <- unique(c(d.2.8))
d.2.8.u[order(d.2.8.u)]
colors <- c(0,d.2.8.u[order(d.2.8.u)])
my_palette <- redgreen(16)

heatmap.3(d.2.8,col=my_palette,Rowv = F,Colv=F,dendrogram="column",
          breaks=colors, density.info="none", trace="none", key=T,symm=T,symkey=F,symbreaks=T,margin=c(9, 13),key.xtickfun=function() {
            breaks <- parent.frame()$breaks
            return(list(
              at=parent.frame()$scale01(c(breaks[1],breaks[length(breaks)])),
              labels=c(as.character(breaks[1]),as.character(breaks[length(breaks)]))
            ))
          })


save.image(file = file.path(output.file.dir,"addMoreBed4JCHeatMap.RData"))
savehistory(file = file.path(output.file.dir,"addMoreBed4JCHeatMap.Rhistory"))

save.image(file = file.path(output.file.dir,"addMoreBed4JCHeatMapSubset.RData"))
savehistory(file = file.path(output.file.dir,"addMoreBed4JCHeatMapSubset.Rhistory"))
