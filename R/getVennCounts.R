##require library graph, RBGL

getVennCounts <- function(..., maxgap = -1L, minoverlap=0L,
                          by=c("region", "feature", "base"), 
                          ignore.strand=TRUE, 
                          connectedPeaks=c("min", "merge", "keepAll")){
  ###check inputs
  PeaksList <- list(...)
  n <- length(PeaksList)
  if(n==1){
    PeaksList <- PeaksList[[1]]
    n <- length(PeaksList)
    names <- names(PeaksList)
    if(is.null(names)) names <- paste("peaks", 1:n, sep="")
  }else{
    ##save dots arguments names
    dots <- substitute(list(...))[-1]
    names <- unlist(sapply(dots, deparse))
  }
  if(n<2){
    stop("Missing required argument Peaks!")
  }
  if(n>5){
    stop("The length of input peaks list should no more than 5")
  }
  connectedPeaks <- match.arg(connectedPeaks)
  by <- match.arg(by)
  if(any(duplicated(names)))
    stop("Same input Peaks detected!")
  venn_cnt <- vennCounts(PeaksList, n=n, names=names, 
                         maxgap=maxgap, minoverlap=minoverlap, by=by,
                         ignore.strand=ignore.strand, 
                         connectedPeaks=connectedPeaks)
  venn_cnt$venn_cnt
}
