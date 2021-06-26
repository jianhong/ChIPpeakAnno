#' Find the overlapped peaks among two or more set of peaks.
#' 
#' Find the overlapping peaks for two or more (less than five) set of peak
#' ranges.
#' 
#' Efficiently perform overlap queries with an interval tree implemented with
#' GRanges.
#' 
#' @aliases findOverlapsOfPeaks overlappingPeaks overlappingPeaks-class
#' @param \dots Objects of \link[GenomicRanges:GRanges-class]{GRanges}: See
#' example below.
#' @param maxgap,minoverlap Used in the internal call to \code{findOverlaps()}
#' to detect overlaps. See
#' \code{?\link[IRanges:findOverlaps-methods]{findOverlaps}} in the
#' \pkg{IRanges} package for a description of these arguments. 
#' If 0 < minoverlap < 1, the function will find overlaps by percentage
#' covered of interval and the filter condition will be set to max covered 
#' percentage of overlapping peaks.
#' @param ignore.strand When set to TRUE, the strand information is ignored in
#' the overlap calculations.
#' @param connectedPeaks If multiple peaks are involved in any group of 
#' connected/overlapping peaks in any input peak list, set it to "merge" will
#' add 1 to the overlapping counts, while set it to "min" will add the minimal 
#' involved peaks in each group of connected/overlapped peaks to the 
#' overlapping counts. Set it to "keepAll" will add the number of involved 
#' peaks for each peak list to the corresponding overlapping counts. 
#' In addition, it will output counts as if connectedPeaks were set to "min".
#' For examples (https://support.bioconductor.org/p/133486/#133603), 
#' if 5 peaks in group1 overlap with 2 peaks in group 2, setting connectedPeaks
#'  to "merge" will add 1 to the overlapping counts; setting it to "keepAll" 
#'  will add 5 peaks to count.group1, 2 to count.group2, and 2 to counts; 
#'  setting it to “min” will add 2 to the overlapping counts.
#' @return return value is An object of overlappingPeaks.  \item{venn_cnt}{an
#' object of VennCounts} \item{peaklist}{a list consists of all overlapping
#' peaks or unique peaks} \item{uniquePeaks}{an object of
#' \link[GenomicRanges:GRanges-class]{GRanges} consists of all unique peaks}
#' \item{mergedPeaks}{an object of \link[GenomicRanges:GRanges-class]{GRanges}
#' consists of all merged overlapping peaks} \item{peaksInMergedPeaks}{an
#' object of \link[GenomicRanges:GRanges-class]{GRanges} consists of all peaks
#' in each samples involved in the overlapping peaks} \item{overlappingPeaks}{a
#' list of data frame consists of the annotation of all the overlapped peaks}
#' \item{all.peaks}{a list of GRanges object which contain the input peaks with
#' formated rownames.}
#' @author Jianhong Ou
#' @seealso \link{annotatePeakInBatch}, \link{makeVennDiagram},
#' \link{getVennCounts}, \link{findOverlappingPeaks}
#' @references 1.Interval tree algorithm from: Cormen, Thomas H.; Leiserson,
#' Charles E.; Rivest, Ronald L.; Stein, Clifford. Introduction to Algorithms,
#' second edition, MIT Press and McGraw-Hill. ISBN 0-262-53196-8
#' 
#' 2.Zhu L.J. et al. (2010) ChIPpeakAnno: a Bioconductor package to annotate
#' ChIP-seq and ChIP-chip data. BMC Bioinformatics 2010,
#' 11:237doi:10.1186/1471-2105-11-237
#' 
#' 3. Zhu L (2013). "Integrative analysis of ChIP-chip and ChIP-seq dataset."
#' In Lee T and Luk ACS (eds.), Tilling Arrays, volume 1067, chapter 4, pp.
#' -19. Humana Press. http://dx.doi.org/10.1007/978-1-62703-607-8_8,
#' http://link.springer.com/protocol/10.1007\%2F978-1-62703-607-8_8
#' @keywords misc
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @examples
#' 
#' peaks1 <- GRanges(seqnames=c(6,6,6,6,5),
#'                  IRanges(start=c(1543200,1557200,1563000,1569800,167889600),
#'                          end=c(1555199,1560599,1565199,1573799,167893599),
#'                          names=c("p1","p2","p3","p4","p5")),
#'                  strand="+")
#' peaks2 <- GRanges(seqnames=c(6,6,6,6,5),
#'                   IRanges(start=c(1549800,1554400,1565000,1569400,167888600),
#'                           end=c(1550599,1560799,1565399,1571199,167888999),
#'                           names=c("f1","f2","f3","f4","f5")),
#'                   strand="+")
#' t1 <- findOverlapsOfPeaks(peaks1, peaks2, maxgap=1000)
#' makeVennDiagram(t1)
#' t1$venn_cnt
#' t1$peaklist
#' t2 <- findOverlapsOfPeaks(peaks1, peaks2, minoverlap = .5)
#' makeVennDiagram(t2)
#' 
#' t3 <- findOverlapsOfPeaks(peaks1, peaks2, minoverlap = .90)
#' makeVennDiagram(t3)
#' 
findOverlapsOfPeaks <- function(..., maxgap=-1L, minoverlap=0L,
                                ignore.strand=TRUE, connectedPeaks=c("keepAll", "min", "merge")){
  ###check inputs
  NAME_conn_string <- "___conn___"
  NAME_short_string <- "__"
  NAME_long_string <- "///"
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
    names <- make.names(unlist(sapply(dots, deparse)))
    names(PeaksList) <- names
  }
  if(any(grepl(NAME_short_string, names))){
    stop(paste("The name of peaks could not contain", NAME_short_string))
  }
  if(n<2){
    stop("Missing required argument Peaks!")
  }
  if(n>5){
    stop("The length of input peaks list should no more than 5")
  }
  connectedPeaks <- match.arg(connectedPeaks)
  if(any(duplicated(names)))
    stop("Same input Peaks detected!")
  PeaksList<-lapply(PeaksList, trimPeakList, by="region", 
                    ignore.strand=ignore.strand, keepMetadata=TRUE)
  venn_cnt <- 
    vennCounts(PeaksList, n=n, names=names, 
               maxgap=maxgap, minoverlap=minoverlap, by="region",
               ignore.strand=ignore.strand, connectedPeaks=connectedPeaks)
  
  outcomes <- venn_cnt$venn_cnt[, 1:n]
  xlist <- do.call(rbind, venn_cnt$xlist)
  xlist <- xlist - 1
  xlist <- xlist[nrow(xlist):1,,drop=FALSE] 
  ## reverse xlist to match the order of names
  xlist <- apply(xlist, 2, base::paste, collapse="")
  if(length(venn_cnt$all)!=length(xlist)) 
    stop("length of 'xlist' and 'all' should be identical.")
  #     all <- do.call(rbind,
  #                    mapply(function(.ele, .id, .gp) cbind(.id, .ele, .gp), 
  #                           venn_cnt$all, 1:length(venn_cnt$all), 
  #                           xlist, SIMPLIFY=FALSE))
  all <- cbind(.id=rep(1:length(venn_cnt$all), sapply(venn_cnt$all, length)), 
               .ele=unlist(venn_cnt$all), 
               .gp=rep(xlist, sapply(venn_cnt$all, length)))
  all.peaks <- venn_cnt$Peaks[all[,2]]
  names(all.peaks) <- gsub(NAME_conn_string,
                           NAME_short_string,
                           names(all.peaks))
  if(!is.null(all.peaks$old_strand_HH)){
    strand(all.peaks) <- all.peaks$old_strand_HH
    all.peaks$old_strand_HH <- NULL
  }
  all.peaks$gpForFindOverlapsOfPeaks <- all[, 1]
  all.peaks$gpType <- all[, 3]
  all.peaks.split <- split(all.peaks, all.peaks$gpType)
  listname <- apply(outcomes, 1, 
                    function(id) 
                      paste(names[as.logical(id)], 
                            collapse=NAME_long_string))
  listcode <- apply(outcomes, 1, base::paste, collapse="")
  listname <- listname[-1]
  listcode <- listcode[-1]
  names(listname) <- listcode
  names(all.peaks.split) <- listname[names(all.peaks.split)]
  peaklist <- lapply(all.peaks.split, reduce, 
                     min.gapwidth=maxgap+1L, with.revmap=TRUE, 
                     ignore.strand=ignore.strand)
  peaklist <- mapply(function(peaks, info){
    revmap <- peaks$revmap
    #         peakNames <- do.call(rbind, 
    #                              mapply(function(id, gp) cbind(id, gp), 
    #                                     revmap, 1:length(revmap), 
    #                                     SIMPLIFY=FALSE))
    peakNames <- cbind(id=unlist(revmap),
                       gp=rep(1:length(revmap), sapply(revmap, length)))
    peaks$peakNames <- CharacterList(split(names(info)[peakNames[, 1]], 
                                           peakNames[, 2]), 
                                     compress=TRUE)
    if(ignore.strand){
      strand <- split(as.character(strand(info))[peakNames[, 1]], 
                      peakNames[, 2])
      strand <- lapply(strand, unique)
      l <- sapply(strand, length)
      strand[l>=2] <- "*"
      strand(peaks) <- unlist(strand)
    }
    peaks$revmap <- NULL
    peaks
  }, peaklist, all.peaks.split, SIMPLIFY=FALSE)
  
  listcode <- strsplit(listcode, "")
  names(listcode) <- listname
  listcode <- sapply(listcode, function(.ele) sum(as.numeric(.ele))==2)
  overlappingPeaks <- sapply(names(listcode)[listcode], function(.ele){
    peakListName <- strsplit(.ele, NAME_long_string, fixed = TRUE)[[1]]
    ps <- PeaksList[peakListName]
    ol <- findOverlaps(query=ps[[1]], subject=ps[[2]], 
                       maxgap=maxgap, minoverlap=minoverlap, 
                       ignore.strand=ignore.strand)
    q <- ps[[1]][queryHits(ol)]
    s <- ps[[2]][subjectHits(ol)]
    cl <- getRelationship(q, s)[,c("insideFeature", "shortestDistance")]
    colnames(cl)[grepl("insideFeature", colnames(cl))] <- "overlapFeature"
    correlation <- cbind(peaks1=names(q), as.data.frame(unname(q)), 
                         peaks2=names(s), as.data.frame(unname(s)), 
                         cl)
    rownames(correlation) <- make.names(paste(names(q), names(s), sep="_"))
    correlation <- correlation[correlation[,"shortestDistance"]<maxgap | 
                                 correlation[,"overlapFeature"] %in% 
                                 c("includeFeature", "inside",
                                   "overlapEnd", "overlapStart",
                                   "overlap"),]
  }, simplify = FALSE)
  PeaksList <- sapply(PeaksList, trimPeakList, by="region",
                      ignore.strand=ignore.strand,
                      keepMetadata=TRUE)
  for(i in 1:n){
    names(PeaksList[[i]]) <- 
      paste(names[i], names(PeaksList[[i]]), sep=NAME_short_string)
  }
  sharedColnames <- Reduce(function(a, b){
    shared <- intersect(colnames(mcols(a)), colnames(mcols(b)))
    if(length(shared)>0){
      shared <- shared[sapply(shared, function(.ele){
        class(mcols(a)[, .ele])==class(mcols(b)[, .ele])
      })]
    }
    if(length(shared)>0) {
      mcols(a) <- mcols(a)[, shared, drop=FALSE]
    }else{
      mcols(a) <- NULL
    }
    a
  }, PeaksList)
  sharedColnames <- colnames(mcols(sharedColnames))
  uniquePeaks <- peaklist[names(peaklist) %in% 
                            listname[!grepl(NAME_long_string, 
                                            listname, fixed = TRUE)]]
  uniquePeaks <- mapply(function(a, b){
    b <- b[unlist(a$peakNames)]
    if(length(sharedColnames)>0){
      mcols(b) <- mcols(b)[, sharedColnames, drop=FALSE]
    }else{
      mcols(b) <- NULL
    }
    b
  }, uniquePeaks, PeaksList[names(uniquePeaks)])
  if(!is(uniquePeaks, "GRangesList"))
    uniquePeaks <- GRangesList(uniquePeaks)
  uniquePeaks <- unlist(uniquePeaks, use.names = FALSE)
  mergedPeaks <- peaklist[names(peaklist) %in% 
                            listname[grepl(NAME_long_string,
                                           listname, fixed = TRUE)]]
  if(!is(mergedPeaks, "GRangesList"))
    mergedPeaks <- GRangesList(mergedPeaks)
  mergedPeaks <- unlist(mergedPeaks, use.names = FALSE)
  peaksInMergedPeaks <- lapply(PeaksList, function(.ele){
    if(length(sharedColnames)>0){
      mcols(.ele) <- mcols(.ele)[, sharedColnames, drop=FALSE]
    }else{
      mcols(.ele) <- NULL
    }
    .ele
  })
  peaksInMergedPeaks <- unlist(GRangesList(peaksInMergedPeaks), 
                               use.names = FALSE)
  peaksInMergedPeaks <- peaksInMergedPeaks[unlist(mergedPeaks$peakNames)]
  structure(list(venn_cnt=venn_cnt$venn_cnt, 
                 peaklist=peaklist, 
                 uniquePeaks=uniquePeaks,
                 mergedPeaks=mergedPeaks,
                 peaksInMergedPeaks=peaksInMergedPeaks,
                 overlappingPeaks=overlappingPeaks,
                 all.peaks=PeaksList), 
            class="overlappingPeaks")
}
