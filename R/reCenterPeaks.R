reCenterPeaks <- function(peaks, width=2000L, ...){
    stopifnot(inherits(peaks, c("annoGR", "GRanges")))
    peaks.center <- start(peaks) + floor(width(peaks)/2)
    peaks.recentered <- peaks
    start(peaks.recentered) <- peaks.center - floor(width/2)
    width(peaks.recentered) <- width
    if(any(start(peaks.recentered)<1)){
      warning("Some start position of the peaks are less than 1!")
    }
    seqLen <- seqlengths(peaks)
    seqLen <- seqLen[!is.na(seqLen)]
    if(length(seqLen)>0){
      peaks.subset <- 
        peaks.recentered[seqnames(peaks.recentered) %in% names(seqLen)]
      if(any(end(peaks.subset)>seqLen[as.character(seqnames(peaks.subset))])){
        warning("Some end position of the peaks are out of bound!")
      }
    }
    peaks.recentered
}