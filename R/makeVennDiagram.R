#' Make Venn Diagram from a list of peaks
#' @description Make Venn Diagram from two or more peak ranges, 
#' Also calculate p-value to  determine whether those peaks 
#' overlap significantly.
#' @param Peaks A list of peaks in \link[GenomicRanges:GRanges-class]{GRanges}
#'  format: See example below.
#' @param NameOfPeaks Character vector to specify the name of Peaks, 
#' e.g., c("TF1", "TF2"). This will be used as label in the Venn Diagram.
#' @param maxgap,minoverlap Used in the internal call to 
#' \code{findOverlaps()} to detect overlaps.
#' See \code{?\link[IRanges:findOverlaps-methods]{findOverlaps}} 
#' in the \pkg{IRanges} package for a description of these arguments.
#' @param totalTest Numeric value to specify the total number of tests 
#' performed to obtain the list of peaks. It should be much larger than 
#' the number of peaks in the largest peak set.
#' @param by "region", "feature" or "base", default = "region". 
#' "feature" means using feature field in the GRanges for calculating overlap, 
#' "region" means using chromosome range for calculating overlap, 
#' and "base" means calculating overlap in nucleotide level.
#' @param ignore.strand Logical: when set to TRUE, the strand information is
#' ignored in the overlap calculations.
#' @param connectedPeaks If multiple peaks involved in overlapping in 
#' several groups, set it to "merge" will count it as only 1, 
#' while set it to "min" will count it as the minimal involved peaks in 
#' any connected peak group. "keepAll" will show all the orginal counts 
#' for each list while the final counts will be same as "min". 
#' "keepFirstListConsistent" will keep the counts consistent with first list.
#' @param method method to be used for p value calculation. 
#' hyperG means hypergeometric test and permutation means \link{peakPermTest}.
#' @param TxDb An object of \link[GenomicFeatures:TxDb-class]{TxDb}.
#' @param plot logical. If TRUE (default), a venn diagram is plotted. 
#' @param \dots Additional arguments to be passed to 
#' \link[VennDiagram:venn.diagram]{venn.diagram}.
#' @details For customized graph options, 
#' please see venn.diagram in VennDiagram package.
#' @return A p.value is calculated by
#' hypergeometric test or permutation test to determine whether the overlaps of 
#' peaks or features are significant.
#' @importFrom VennDiagram venn.diagram
#' @importFrom grid textGrob gList grid.newpage
#' @export
#' @author Lihua Julie Zhu, Jianhong Ou
#' @seealso \link{findOverlapsOfPeaks}, 
#' \link[VennDiagram:venn.diagram]{venn.diagram}, \link{peakPermTest}
#' @examples 
#' if (interactive()){
#' peaks1 <- GRanges(seqnames=c("1", "2", "3"),
#'                   IRanges(start=c(967654, 2010897, 2496704),
#'                           end=c(967754, 2010997, 2496804), 
#'                           names=c("Site1", "Site2", "Site3")),
#'                   strand="+",
#'                   feature=c("a","b","f"))
#' peaks2 = GRanges(seqnames=c("1", "2", "3", "1", "2"), 
#'                  IRanges(start = c(967659, 2010898,2496700,
#'                                    3075866,3123260),
#'                          end = c(967869, 2011108, 2496920, 
#'                                  3076166, 3123470),
#'                          names = c("t1", "t2", "t3", "t4", "t5")), 
#'                  strand = c("+", "+", "-", "-", "+"), 
#'                  feature=c("a","b","c","d","a"))
#' makeVennDiagram(list(peaks1, peaks2), NameOfPeaks=c("TF1", "TF2"),
#'                 totalTest=100,scaled=FALSE, euler.d=FALSE, 
#'                 fill=c("#009E73", "#F0E442"), # circle fill color
#'                 col=c("#D55E00", "#0072B2"), #circle border color
#'                 cat.col=c("#D55E00", "#0072B2"))
#' 
#' makeVennDiagram(list(peaks1, peaks2), NameOfPeaks=c("TF1", "TF2"),
#'                 totalTest=100, 
#'                 fill=c("#009E73", "#F0E442"), # circle fill color
#'                 col=c("#D55E00", "#0072B2"), #circle border color
#'                 cat.col=c("#D55E00", "#0072B2"))
#' 
#' ###### 4-way diagram using annotated feature instead of chromosome ranges
#' 
#' makeVennDiagram(list(peaks1, peaks2, peaks1, peaks2), 
#'                 NameOfPeaks=c("TF1", "TF2","TF3", "TF4"), 
#'                 totalTest=100, by="feature",
#'                 main = "Venn Diagram for 4 peak lists",
#'                 fill=c(1,2,3,4))
#' }
#' @keywords graph
#' 

makeVennDiagram <- function(Peaks, NameOfPeaks, maxgap=-1L, minoverlap=0L,
                            totalTest, by=c("region", "feature", "base"), 
                            ignore.strand=TRUE, 
                            connectedPeaks=c("min", "merge", "keepAll", 
                                             "keepFirstListConsistent"), 
                            method=c("hyperG", "permutation"), 
                            TxDb, plot=TRUE,
                            ...){
  ###Functions to be used
  getCountsList<-function(counts){
    CountsList <- list()
    cnt <- 1
    for(i in 1:length(counts)){
      CountsList[[i]] <- seq(cnt, length.out=counts[i])
      cnt <- cnt+counts[i]
    }
    CountsList
  }
  getVennList<-function(a, NameOfPeaks, CountsList){
    .x<-lapply(NameOfPeaks, function(.ele, cl, a){
      .y<-c()
      for(i in 1:nrow(a)){
        if(a[i,.ele]!=0) .y <- c(.y, cl[[i]])
      }
      .y
    }, CountsList, a)
    names(.x)<-NameOfPeaks
    .x
  }
  plotVenn<-function(venn_cnt, vennx, otherCounts=NULL, ...){
    dots <- list(...)
    if(is.null(dots$cat.cex)) dots$cat.cex = 1
    if(is.null(dots$cat.col)) dots$cat.col = "black"
    if(is.null(dots$cat.fontface)) dots$cat.fontface = "plain"
    if(is.null(dots$cat.fontfamily)) dots$cat.fontfamily = "serif"
    
    op=par(mar=c(0,0,0,0))
    on.exit(par(op))
    grid.newpage()
    dots$x <- vennx
    dots <- c(dots, filename=list(NULL))
    venngrid <- do.call(venn.diagram, dots)
    unlink(dir(pattern="^VennDiagram[0-9_\\-]+.log$")) ## delete the log file
    if(grepl("^count\\.", colnames(venn_cnt)[ncol(venn_cnt)]) && 
           connectedPeaks=="keepAll"){
        n <- which(colnames(venn_cnt)=="Counts")-1
        venn_cnt.Counts.gt.1 <- rowSums(venn_cnt[ ,1:n]) > 1
        counts <- venn_cnt[-1,"Counts"]
        counts <- counts[counts!=0]
        if(!any(duplicated(counts))){## there is no duplicated counts number
            w <- 1
            for(i in 1:length(venngrid)){
                if(inherits(venngrid[[i]], what="text") && w<nrow(venn_cnt)){
                    w <- w+1
                    cnt <- as.numeric(venngrid[[i]]$label)
                    if(!is.na(cnt) && cnt!=0){
                        j <- which(venn_cnt[,n+1]==cnt)
                        if(venn_cnt.Counts.gt.1[j]){
                            labels <- venn_cnt[j, -(1:(n+1))]
                            if(nrow(venn_cnt)==4 && sum(venn_cnt[,4]) < 
                                   sum(venn_cnt[,5])){
                                labels <- rev(labels)
                            }
                            labels <- paste(as.character(labels[labels>0]),
                                            collapse="/", sep="")
                            if(labels!="") venngrid[[i]]$label <- 
                                paste(venngrid[[i]]$label, "(", 
                                      labels, ")", sep="")
                        }
                    }
                }
            }
        }else{##otherwise guess the order of the number, 
              ##Bugs, the map order maybe wrong
            j <- 1
            map <- switch(as.character(nrow(venn_cnt)),
                          '4'=c(3, 4, 2),
                          '8'=c(5, 7, 3, 6, 8, 4, 2),
                          '16'=c(3, 4, 2, 11, 12, 16, 8, 6, 9, 10, 
                                 14, 15, 7, 5, 13),
                          '32'=c(17, 9, 5, 3, 2, 6, 18, 19, 25, 10, 13, 21, 7, 
                                 11, 4, 8, 22, 20, 27, 26, 30, 29, 23, 
                                 15, 12, 16, 24, 28, 30, 31, 32))
            for(i in 1:length(venngrid)){
                if(inherits(venngrid[[i]], what="text") && j<length(map)){
                    while(j < length(map) && as.numeric(venngrid[[i]]$label) != 
                              venn_cnt[map[j],n+1]) j <- j+1
                    if(as.numeric(venngrid[[i]]$label) == venn_cnt[map[j],
                                                                   n+1] && 
                           j <= length(map)){
                        if(venn_cnt.Counts.gt.1[map[j]]){
                            labels <- venn_cnt[map[j], -(1:(n+1))]
                            labels <- paste(as.character(labels[labels>0]), 
                                            collapse="/", sep="")
                            if(labels!="") 
                                venngrid[[i]]$label <- 
                                paste(venngrid[[i]]$label, "(", 
                                      labels, ")", sep="")
                        }
                        j <- j+1
                    }
                }
            }
        }
    }
    if(!is.null(otherCounts)){
      tmp <- textGrob(label=paste("others:", otherCounts), x=0.9, y=0.1, 
                      gp=gpar(col = "black", cex = dots$cat.cex, 
                              fontface = dots$cat.fontface, 
                              fontfamily = dots$cat.fontfamily))
      venngrid <- gList(venngrid, tmp)
    }
    grid.draw(venngrid)
  }
  getPval<-function(venn_cnt, totalTest){
    n <- which(colnames(venn_cnt)=="Counts") - 1 ##ncol(venn_cnt)-1
    s <- apply(venn_cnt[,1:n], 1, sum)
    venn_cnt_s <- venn_cnt[s==2, , drop=FALSE]
	  for(i in 1:nrow(venn_cnt_s)){
		  venn_cnt_s[i,n+1] <- 
              sum(venn_cnt[as.numeric(venn_cnt[,1:n]%*%venn_cnt_s[i,1:n])==2,
                           n+1,drop=TRUE])
	  }
    cnt <- venn_cnt[,n+1]
    cnt_s <- apply(venn_cnt[,1:n], 2, function(.ele){
      sum(cnt[as.logical(.ele)])
    })
    p.value <- apply(venn_cnt_s, 1, function(.ele){
      ab <- cnt_s[as.logical(.ele[1:n])]
      a <- ab[1]
      b <- ab[2]
      a.and.b <- .ele[(n+1)]
      phyper <- phyper(a.and.b-1, b, totalTest-b, a, 
                       lower.tail=FALSE, log.p=FALSE)
    })
    cbind(venn_cnt_s[, 1:n, drop=FALSE], pval=p.value)
  }
  ###check inputs
  method <- match.arg(method)
  if (missing(totalTest) && method=="hyperG")
  {
    message("Missing totalTest! totalTest is required for HyperG test. 
If totalTest is missing, pvalue will be calculated by estimating 
the total binding sites of encoding region of human.
totalTest = humanGenomeSize * (2%(codingDNA) + 
             1%(regulationRegion)) / ( 2 * averagePeakWidth )
          = 3.3e+9 * 0.03 / ( 2 * averagePeakWidth)
          = 5e+7 /averagePeakWidth")
  }
  if (missing(Peaks))
  {
    stop("Missing Peaks which is a list of peaks in GRanges!")
  }
  connectedPeaks <- match.arg(connectedPeaks)
  stopifnot(is.logical(plot))
  by <- match.arg(by)
  n1 = length(Peaks)
  olout_flag <- FALSE
  if(is(Peaks, "overlappingPeaks")){
    if(!is.null(Peaks$venn_cnt) && is(Peaks$venn_cnt, "VennCounts")){
        venn_cnt <- Peaks$venn_cnt
        n1 <- which(colnames(venn_cnt)=="Counts") - 1 ##ncol(venn_cnt)-1
        if(missing(NameOfPeaks)){
            NameOfPeaks <- colnames(venn_cnt)[1:n1]
        }
        olout_flag <- TRUE
    }else{
        stop("Input is object of overlappingPeaks. 
             But it does not have venn counts data.")
    }
  }
  if (missing(NameOfPeaks) ||  mode(NameOfPeaks) != "character")
  {
    warning("Missing required character vector NameOfPeaks. 
            NameOfPeaks will be extract from the names of input Peaks.")
    NameOfPeaks <- names(Peaks)
    if(is.null(NameOfPeaks)){
      NameOfPeaks <- paste("peaks", 1:length(Peaks), sep="")
    }
  }
  n2 = length(NameOfPeaks)
  if (n1 <2)
  {
    stop("The number of element in Peaks is less than 1, 
         need at least 2 peak ranges!")
  }
  if (n1 > 5) stop("The length of input peaks list should no more than 5")
  if (n1 > n2)
  {
     stop("The number of element in NameOfPeaks is less 
          than the number of elements in Peaks, need to be equal!")
  }
  if (n1 < n2)
  {
    warning("The number of element in NameOfPeaks is larger than 
            the number of elements in Peaks! 
            NameOfPeaks will be cut by the length of Peaks")
    NameOfPeaks <- NameOfPeaks[1:n1]
  }
  NameOfPeaks <- make.names(NameOfPeaks, unique=TRUE, allow_=TRUE)
  if(!olout_flag){
      for(i in seq_len(n1)){
        if (!inherits(Peaks[[i]], c("GRanges")))
        {
          err.msg = paste("Element", i, 
                          "in Peaks is not a valid GRanges object", sep=" ")
          stop(err.msg)
        }
        if (is.null(names(Peaks[[i]]))) names(Peaks[[i]]) <-  
            formatC(1:length(Peaks[[i]]), 
                    width=nchar(length(Peaks[[i]])), 
                    flag="0")
      }
      if(!missing(totalTest)){
          if(totalTest < max(sapply(Peaks, length))){
              stop("totaltest specifies the total number of possible peaks 
                   in the testing space. It should be larger than the largest
                   peak numbers in the input sets. Please see more details 
                   at http://pgfe.umassmed.edu/ChIPpeakAnno/FAQ.html")
          }
      }
      names(Peaks) <- NameOfPeaks
      venn_cnt <- 
        getVennCounts(Peaks, maxgap=maxgap, 
                      minoverlap=minoverlap, 
                      by=by, ignore.strand=ignore.strand, 
                      connectedPeaks=
                        ifelse(connectedPeaks=="keepFirstListConsistent", 
                               "keepAll", connectedPeaks))
  }
  colnames(venn_cnt)[1:n1] <- NameOfPeaks
  venn_cnt1 <- venn_cnt
  if(connectedPeaks=="keepFirstListConsistent"){ ## keepAll to getVennCounts
    if(!grepl("^count\\.", colnames(venn_cnt)[ncol(venn_cnt)])){
      warning("If connectedPeaks set to keepFirstListConsistent,",
              "connectedPeaks of findOverlapsOfPeaks must be keepAll.",
              "Setting connectedPeaks to default.")
    }else{
      venn_cnt1[venn_cnt[, 1]==1, "Counts"] <- venn_cnt[venn_cnt[, 1]==1, n1+2]
    }
    connectedPeaks <- "min"
  }
  Counts <- getCountsList(venn_cnt1[,"Counts"])
  vennx <- getVennList(venn_cnt1, NameOfPeaks, Counts)
  if(method=="hyperG"){
      if(!missing(totalTest)){
          otherCount <- totalTest - sum(venn_cnt[, "Counts"])
          venn_cnt[1, "Counts"] <- otherCount
          p.value <- getPval(venn_cnt, totalTest)
      }else{
          otherCount <- venn_cnt[1, "Counts"]
          if(is(Peaks, "overlappingPeaks")){
              averagePeakWidth = median(unlist(lapply(Peaks$peaklist, width)))
          }else{
              averagePeakWidth = median(unlist(lapply(Peaks, width))) 
          }
          
          p.value <- getPval(venn_cnt, round(5e+7/averagePeakWidth))
      }
  }else{
      ##method=="permutation"
      p.value <- NA
      if(olout_flag){
          warning("Input is an object of overlappingPeaks.
                  Can not do permutation test. Please try ?peakPermTest")
      }else{
          if(missing(TxDb)){
              warning("TxDb is missing. Please try ?peakPermTest later.")
          }else{
              n <- which(colnames(venn_cnt)=="Counts") - 1 ##ncol(venn_cnt)-1
              s <- apply(venn_cnt[,1:n], 1, sum)
              venn_cnt_s <- venn_cnt[s==2, , drop=FALSE]
              p.value <- apply(venn_cnt_s, 1, function(.ele){
                  AB <- Peaks[as.logical(.ele)]
                  pt <- peakPermTest(AB[[1]], AB[[2]], TxDb=TxDb, ...)
                  pt$pval
              })
              p.value <- cbind(venn_cnt_s[, 1:n, drop=FALSE], pval=p.value)
          }
      }
      if(!missing(totalTest)){
          otherCount <- totalTest - sum(venn_cnt[, "Counts"])
          venn_cnt[1, "Counts"] <- otherCount
      }else{
          otherCount <- venn_cnt[1, "Counts"]
      }
  }
  if(otherCount==0) otherCount <- NULL
  if(plot) {
    plotVenn(venn_cnt1, vennx, otherCounts=otherCount, ...)
  }
  return(list(p.value=p.value, vennCounts=venn_cnt))
}
