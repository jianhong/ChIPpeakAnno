##require library graph, RBGL

#' Obtain Venn Counts for Venn Diagram, internal function for makeVennDigram
#' 
#' Obtain Venn Counts for peak ranges using chromosome ranges or feature field,
#' internal function for makeVennDigram
#' 
#' 
#' @param \dots Objects of \link[GenomicRanges:GRanges-class]{GRanges}. See
#' example below.
#' @param maxgap,minoverlap Used in the internal call to \code{findOverlaps()}
#' to detect overlaps. See
#' \code{?\link[IRanges:findOverlaps-methods]{findOverlaps}} in the
#' \pkg{IRanges} package for a description of these arguments.
#' @param by region, feature or base, default region. feature means using
#' feature field in the GRanges for calculating overlap, region means using
#' chromosome range for calculating overlap, and base means using calculating
#' overlap in nucleotide level.
#' @param ignore.strand When set to TRUE, the strand information is ignored in
#' the overlap calculations.
#' @param connectedPeaks If multiple peaks involved in overlapping in several
#' groups, set it to "merge" will count it as only 1, while set it to "min"
#' will count it as the minimal involved peaks in any concered groups
#' @return \item{vennCounts}{vennCounts objects containing counts for Venn
#' Diagram generation, see details in limma package vennCounts}
#' @author Jianhong Ou
#' @seealso \link{makeVennDiagram}, \link{findOverlappingPeaks}
#' @keywords misc
#' @export
#' @examples
#' 
#' if(interactive() || Sys.getenv("USER")=="jianhongou"){
#' peaks1 = GRanges(seqnames=c("1", "2", "3"), 
#'                  IRanges(start = c(967654, 2010897, 2496704), 
#'                             end = c(967754, 2010997, 2496804), 
#'                             names = c("Site1", "Site2", "Site3")),  
#'                    strand=as.integer(1), 
#'                    feature=c("a","b", "c"))
#'   peaks2 = 
#'       GRanges(seqnames= c("1", "2", "3", "1", "2"), 
#'                     IRanges(start=c(967659, 2010898, 2496700, 3075866, 3123260), 
#'                          end=c(967869, 2011108, 2496920, 3076166, 3123470), 
#'                          names = c("t1", "t2", "t3", "t4", "t5")), 
#'                     strand = c(1L, 1L, -1L,-1L,1L), 
#'                     feature=c("a","c","d","e", "a"))
#'     getVennCounts(peaks1,peaks2)
#'     getVennCounts(peaks1,peaks2, by="feature")
#'     getVennCounts(peaks1, peaks2, by="base")
#' }
#' 
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
