TEAD.bed <- "/Users/aiminyan/Aimin/DropboxUmass/Aimin/Project/Dimpi_Mercurio_lab/ENCFF216XPT.bed"

bed.files <- c(TEAD.bed)

bed.in<-lapply(1:length(bed.files),function(u,bed.files){
  
  x <- bed.files[[u]]
  
  if(u == 1){peaks=read.table(x,header=F)}else{peaks=read.table(x,header=F,skip=1)}
  
  colnames(peaks)[1:3]= c("chr","start","end")
  peaks=toGRanges(peaks)
  peaks
},bed.files)

bed.4.gene.density.1.anno <- getAnnotatedGene(bed.in,"Hs")

grep("TRPC6",bed.4.gene.density.1.anno[[1]]$gene_name)

grep("ENSG00000137672",bed.4.gene.density.1.anno[[1]]$feature)
     
grep("SRD5A1P1",bed.4.gene.density.1.anno[[1]]$gene_name)

bed.4.gene.density.1.anno[[1]][2586,]

which(start(bed.4.gene.density.1.anno[[1]])==101337181,)

grep("DYNC2H1",bed.4.gene.density.1.anno[[1]]$gene_name)

grep("YAP1",bed.4.gene.density.1.anno[[1]]$gene_name)
bed.4.gene.density.1.anno[[1]][903,]

grep("LUC7L",bed.4.gene.density.1.anno[[1]]$gene_name)

bed.4.gene.density.1.anno[[1]][903,]

library(EnsDb.Hsapiens.v75)
annoData <- annoGR(EnsDb.Hsapiens.v75)
annotatePeak = annotatePeakInBatch(bed.in[[1]], AnnotationData=annoData)

# "ENSG00000137672": TRPC6 
grep("ENSG00000137672",annotatePeak$feature)

# "ENSG00000137693": YAP1
grep("ENSG00000137693",annotatePeak$feature)

annotatePeak[c(3677,3678),]

# "overlapping" to TRPC6("ENSG00000137672")
AnnotatePeak.overlapping = annotatePeakInBatch(bed.in[[1]], AnnotationData=annoData,output = "overlapping")
grep("ENSG00000137672",AnnotatePeak.overlapping$feature)
AnnotatePeak.overlapping[c(4020,4021,4022),]

# does not work
AnnotatePeak.nearestBiDirectionalPromoters = annotatePeakInBatch(bed.in[[1]], AnnotationData=annoData,output = "nearestBiDirectionalPromoters")

# "nearestBiDirectionalPromoters" to TRPC6("ENSG00000137672")
AnnotatePeak.nearestBiDirectionalPromoters = annoPeaks(bed.in[[1]],annoData,output = "nearestBiDirectionalPromoters")
grep("ENSG00000137672",AnnotatePeak.nearestBiDirectionalPromoters$feature)

output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/Aimin/Project/Dimpi_Mercurio_lab"

goodGR <- trim(AnnotatePeak.overlapping[c(4020,4021,4022),])
overlaps.trimmed<-goodGR

library(BSgenome.Hsapiens.UCSC.hg19)

seq<-getAllPeakSequence(overlaps.trimmed,genome=Hsapiens,upstream=0, downstream=0)
write2FASTA(seq, file.path(output.file.dir,"TRPC6.fa"))

save.image(file = file.path(output.file.dir,"Dimpi_Mercurio_lab.RData"))
savehistory(file = file.path(output.file.dir,"Dimpi_Mercurio_lab.Rhistory"))



