files <- getSampleFiles()
print(files)

peak <- readPeakFile(files[[4]])
peak

data("tagMatrixList")
tagMatrix <- tagMatrixList[[4]]

tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")

peakAnno <- annotatePeak(files[[4]], tssRegion=c(-3000, 3000), 
                         TxDb=txdb, annoDb="org.Hs.eg.db")