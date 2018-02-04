#' Plot the cumulative percentage tag allocation in sample
#' 
#' Plot the difference between the cumulative percentage tag allocation in 
#' paired samples.
#' 
#' @param bamfiles Bam file names.
#' @param gr An object of \link[GenomicRanges:GRanges-class]{GRanges}
#' @param input Which file name is input. default 1.
#' @param binWidth The width of each bin.
#' @param ... parameter for \link[GenomicAlignments]{summarizeOverlaps}.
#' @import SummarizedExperiment
#' @import GenomicAlignments
#' @import GenomicRanges
#' @import S4Vectors
#' @importFrom BiocGenerics which.min, which.max
#' @importFrom graphics abline axis legend matlines par plot
#' @export
#' @return A list of data.frame with the cumulative percentages.
#' @references Normalization, bias correction, and peak calling for ChIP-seq
#' Aaron Diaz, Kiyoub Park, Daniel A. Lim, Jun S. Song
#' Stat Appl Genet Mol Biol. Author manuscript; 
#' available in PMC 2012 May 3.Published in final edited form as: 
#' Stat Appl Genet Mol Biol. 2012 Mar 31; 11(3): 10.1515/1544-6115.1750
#'  /j/sagmb.2012.11.issue-3/1544-6115.1750/1544-6115.1750.xml. 
#'  Published online 2012 Mar 31.  doi: 10.1515/1544-6115.1750
#'  PMCID: PMC3342857
#' @examples
#' \dontrun{
#' path <- system.file("extdata", "reads", package="MMDiffBamSubset")
#' files <- dir(path, "bam$", full.names = TRUE)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' gr <- as(seqinfo(Hsapiens)["chr1"], "GRanges")
#' cumulativePercentage(files, gr)
#' } 
#'  

cumulativePercentage <- function(bamfiles, gr, input=1, binWidth=1e3,
                                 ...){
  stopifnot(is(gr, "GRanges"))
  tileTargetRegions <- tileGRanges(gr, windowSize = binWidth, step = binWidth)
  se <- summarizeOverlaps(features=tileTargetRegions, reads=bamfiles, ...)
  ## resample
  sampleName <- colnames(assays(se)[[1]])[-1*input]
  input <- colnames(assays(se)[[1]])[input]
  sigBin <- lapply(sampleName, function(.n){
    sig <- assays(se)[[1]][, c(input, .n)]
    sig[order(sig[, .n], sig[, input]), ]
  })
  sigCumsum <- lapply(sigBin, function(.ele){
    .ele <- apply(.ele, 2, cumsum)
    .ele <- cbind(Rank=seq.int(nrow(.ele)), .ele)
    sweep(.ele, MARGIN = 2, STATS = .ele[nrow(.ele),], FUN = `/`)
  })
  pin <- par("pin")
  if(pin[2]>0){
    ratio <- 2^round(diff(log2(pin)))
    n <- length(sampleName)
    ncol <- ceiling(sqrt(n/ratio))
    nrow <- ceiling(n/ncol)
    op <- par(mfrow=c(nrow, ncol), pty="s")
    on.exit(par(op))
    for(i in 1:n){
      ## plot
      plot(c(0, 1), c(0, 1), type="n", 
           xlab="% of bins", ylab="% of tags", 
           main=sampleName[i])
      matlines(x=sigCumsum[[i]][, 1], 
               y=sigCumsum[[i]][, -1],
               lty = 1)
      zero <- which(sigCumsum[[i]][, 2]>1/binWidth)
      maxDist <- which.max(abs(sigCumsum[[i]][, 2] - sigCumsum[[i]][, 3]))
      if(length(zero)>0 | length(maxDist)>0){
        x <- numeric(0)
        if(length(zero)>0) x <- c(x, zero[1])
        if(length(maxDist)>0) x <- c(x, maxDist[1])
        x.tick <- sigCumsum[[i]][x, 1]
        abline(v=x.tick, col="yellowgreen", lty=3)
        axis(3, at=x.tick, labels = formatC(x.tick, digits=2), 
             lwd = 0, lwd.ticks = 1)
      }
      legend("topleft", legend = colnames(sigCumsum[[i]])[-1], 
             col=1:6, lty = 1, pch = NA, box.col = NA)
    }
  }
  return(invisible(sigCumsum))
}