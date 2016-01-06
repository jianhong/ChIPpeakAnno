reCenterPeaks <- function(peaks, width=2000L, ...){
    stopifnot(inherits(peaks, c("annoGR", "GRanges")))
    peaks.center <- start(peaks) + floor(width(peaks)/2)
    peaks.recentered <- peaks
    start(peaks.recentered) <- peaks.center - floor(width/2)
    width(peaks.recentered) <- width
    peaks.recentered
}