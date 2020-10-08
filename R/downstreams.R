#' Get downstream coordinates
#' 
#' Returns an object of the same type and length as x containing downstream 
#' ranges. The output range is defined as
#' 
#' (end(x) - upstream) to (end(x) + downstream -1) 
#' 
#' for ranges on the + and * 
#' strand, and as 
#' 
#' (start(x) - downstream + 1) to (start(x) + downstream)
#' 
#' for ranges on the - strand.
#' 
#' Note that the returned object might contain out-of-bound ranges.
#' 
#' @param gr A GenomicRanges object
#' @param upstream,downstream non-negative interges.
#' @return A GenomicRanges object
#' @export
#' @import GenomicRanges
#' @examples 
#' gr <- GRanges("chr1", IRanges(rep(10, 3), width=6), c("+", "-", "*"))
#' downstreams(gr, 2, 2)
#'
downstreams <- function(gr, upstream, downstream){
  stopifnot(is(gr, "GRanges"))
  stopifnot(identical(levels(strand(gr)), c("+", "-", "*")))
  gr.rev <- gr
  levels(strand(gr.rev)) <- c("-", "+", "*")
  out <- promoters(gr.rev, upstream = downstream, downstream = upstream)
  levels(strand(out)) <- c("-", "+", "*")
  stopifnot(identical(levels(strand(out)), c("+", "-", "*")))##will auto Correct
  out
}
