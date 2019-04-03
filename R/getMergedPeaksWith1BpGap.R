chr.bed <- bed.in[[6]]

getRestedRegionsFromGR <- function(NAD.peaks.nonXL.in) {
  chr.bed <- NAD.peaks.nonXL.in
  chr.bed.temp <- chr.bed
  start(chr.bed.temp) <- min(start(chr.bed))
  end(chr.bed.temp) <- max(end(chr.bed))
  chr.bed.temp.1 <- unique(chr.bed.temp)
  NAD.peaks.nonXL.rest <- GenomicRanges::(chr.bed.temp.1,chr.bed)
  NAD.peaks.nonXL.rest
}

bed.in.6.X <- bed.in[[6]][which(seqnames(bed.in[[6]])=="chr19")]

XX.3 <- bed.in.3.X[which(start(bed.in.3.X)==3050001)]
YY.3 <- bed.in.3.X[which(start(bed.in.3.X)==3100001)]

XX.YY.3 <- c(XX.3,YY.3)


getMergedPeaksWith1BpGap <- function(XX.YY.3) {
  hits <- GenomicRanges::findOverlaps(XX.YY.3,maxgap=1L, minoverlap=0L)
  net.hits <- hits[!(isSelfHit(hits) | isRedundantHit(hits))]
  XX.YY.3.temp <- XX.YY.3
  start(XX.YY.3.temp) <- start((XX.YY.3[queryHits(net.hits)]))
  end(XX.YY.3.temp) <- end((XX.YY.3[subjectHits(net.hits)]))                           
  XX.YY.3.temp.merged <- unique(GenomicRanges::sort(XX.YY.3.temp))
  XX.YY.3.temp.merged
}

peaks1 <- GRanges(seqnames=c("1", "1", "1","1"),
                  IRanges(start=c(3,8, 15, 60),
                          end=c(7,14, 50, 100)))

peaks2 <- GRanges(seqnames=c("1"),
                  IRanges(start=min(start(peaks1)),
                          end=max(end(peaks1))))

GenomicRanges::setdiff(peaks2,peaks1,fill.gap=T)

hits <- GenomicRanges::findOverlaps(peaks1,maxgap=3, minoverlap=0)
net.hits <- hits[!(isSelfHit(hits) | isRedundantHit(hits))]


XX.YY.3
getMergedPeaksWith1BpGap(XX.YY.3)
getContinuousRegion(XX.YY.3,maxGap=1,minGap=1)
getMergedPeaksWith1BpGap(bed.in.6.X)

bed.in.6.X.sorted <- unique(GenomicRanges::sort(bed.in.6.X))

non.olp.index <- getNonOverLappingPeakIndex(bed.in[[6]])
bed.in.6.non.overlap <- bed.in[[6]][non.olp.index]
bed.in.6.overlap <- bed.in[[6]][-non.olp.index]

bed.in.6.non.overlap.2 <- getMergedPeaksWith1BpGap(bed.in.6.X)

hits <- GenomicRanges::findOverlaps(XX.YY.3,maxgap=1L, minoverlap=0L)
net.hits <- hits[!(isSelfHit(hits) | isRedundantHit(hits))]
XX.YY.3.temp <- XX.YY.3

GenomicRanges::reduce(bed.in.6.X)

hits <- GenomicRanges::findOverlaps(bed.in.6.non.overlap,maxgap=1L, minoverlap=0L)
net.hits <- hits[!(isSelfHit(hits) | isRedundantHit(hits))]

bed.in[[6]][queryHits(net.hits)]

bed.in.6.X.merged <- getMergedPeaksWith1BpGap(bed.in.6.X.sorted) 


bed.in.3.X <- bed.in[[3]][which(seqnames(bed.in[[3]])=="chr19")]
res.Hidden_Domains <- getRestedRegionsFromGR(bed.in[[6]])

bed.in[8]
GenomicRanges::reduce(bed.in[[8]])

bed.in[[4]]
GenomicRanges::reduce(bed.in[[4]])

bed.in[[5]]
GenomicRanges::reduce(bed.in[[5]])

null <- lapply(1:length(bed.in), function(u,bed.in){
  cat(names(bed.in)[u],"\n")
  cat(length(bed.in[[u]]),"\n")
  cat(length(GenomicRanges::reduce(bed.in[[u]])),"\n")
},bed.in)


bed.in.refined <- lapply(1:length(bed.in), function(u,bed.in){
  XX <- GenomicRanges::reduce(bed.in[[u]])
  XX
},bed.in)

names(bed.in.refined) <- names(bed.in)

output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/NADSplittingRegion/RefinedPeaks"

null <- lapply(1:length(bed.in.refined), function(u,bed.in.refined){
  
  output.file.name <- paste0(names(bed.in.refined)[u],"_refined")
  
  cat(output.file.name,"\n")
  orderPeakAndOutPut(as.data.frame(bed.in.refined[[u]]),output.file.dir,output.file.name,outHeader= FALSE)
  
},bed.in.refined)

bed.in.refined[[6]][which(start(bed.in.refined[[6]])==47750000)]

bed.in[[6]][which(start(bed.in[[6]])==47750000)]


