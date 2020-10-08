#' Genomic Element data for upset plot
#' 
#' Prepare data for upset plot for genomic element distribution
#' 
#' @details The data will be calculated by for each breaks. 
#' No precedence will be considered. 
#' @param peaks peak list, \link[GenomicRanges:GRanges-class]{GRanges} object or
#' a \link[GenomicRanges:GRangesList-class]{GRangesList}.
#' @param TxDb an object of \code{\link[GenomicFeatures:TxDb-class]{TxDb}}
#' @param seqlev sequence level should be involved. 
#' Default is all the sequence levels in intersect of peaks and TxDb.
#' @param ignore.strand logical. Whether the strand of the input ranges
#' should be ignored or not. Default=TRUE
#' @param breaks list. A list for labels and sets for the genomic elements.
#' The element could be an S4 method for signature 'TxDb' or a numeric vector
#' with length of 4. The three numbers are 
#' c(upstream point, downstream point, promoter (-1) or downstream (1),
#'  remove gene body or not (1: remove, 0: keep)).
#' @export
#' @importFrom GenomicFeatures intronsByTranscript exons fiveUTRsByTranscript 
#' threeUTRsByTranscript genes cds promoters
#' @return list of data for plot.
#' @examples 
#' if (interactive() || Sys.getenv("USER")=="jianhongou"){
#'   data(myPeakList)
#'   if(require(TxDb.Hsapiens.UCSC.hg19.knownGene)){
#'   seqinfo(myPeakList) <- 
#'   seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene)[seqlevels(myPeakList)]
#'   myPeakList <- GenomicRanges::trim(myPeakList)
#'   myPeakList <- myPeakList[width(myPeakList)>0]
#'   x <- genomicElementUpSetR(myPeakList, 
#'     TxDb.Hsapiens.UCSC.hg19.knownGene)
#'   library(UpSetR)
#'   upset(x$plotData, nsets=13, nintersects=NA)
#'   }
#' }
genomicElementUpSetR <- 
  function(peaks, TxDb, seqlev, ignore.strand=TRUE,
           breaks = 
             list("distal_upstream"=c(-100000, -10000, -1, 1), 
                  "proximal_upstream"=c(-10000,-5000, -1, 1),
                  "distal_promoter"=c(-5000, -2000, -1, 1),
                  "proximal_promoter"=c(-2000, 200, -1, 0),
                  "5'UTR"=fiveUTRsByTranscript,
                  "3'UTR"=threeUTRsByTranscript,
                  "CDS"=cds,
                  "exon"=exons,
                  "intron"=intronsByTranscript,
                  "gene_body"=genes,
                  "immediate_downstream"=c(0, 2000, 1, 1),
                  "proximal_downstream"=c(2000, 5000, 1, 1),
                  "distal_downstream"=c(5000, 100000, 1, 1))){
    stopifnot("peaks must be an object of GRanges or GRangesList"=
                inherits(peaks, c("GRanges", "GRangesList")))
    if(is(peaks, "GRanges")){
      n <- deparse(substitute(peaks))
      peaks <- GRangesList(peaks)
      names(peaks) <- n
      isGRanges <- TRUE
    }else{
      isGRanges <- FALSE
    }
    stopifnot("TxDb must be an object of TxDb"=is(TxDb, "TxDb"))
    
    defaultW <- getOption("warn")
    options(warn = -1)
    on.exit({warn = defaultW})
    
    seql <- seqlevelsStyle(peaks)
    
    ## set annotation
    suppressMessages(g <- genes(TxDb, single.strand.genes.only=TRUE))
    
    anno <- lapply(breaks, function(.ele){
      stopifnot("Elements of breaks must be function or numeric(4)"=
                  is.function(.ele)||is.numeric(.ele))
      if(is.numeric(.ele)){
        stopifnot("Elements of breaks must be function or numeric(4)"=
                    length(.ele)==4)
        ## upstream or downstream
        ups_dws <- function(n){
          ifelse(n[3]==-1,
                 ifelse(all(n[1:2]<0), "uu", 
                        ifelse(all(n[1:2]>0), "pp", "ups")),
                 ifelse(all(n[1:2]>0), "dd", 
                        ifelse(all(n[1:2]<0), "ww", "dws")))
        }
        x <- ups_dws(.ele)
        fil <- .ele[4]
        .ele <- .ele[1:2]
        .ele <- switch(x,
                      "uu"={
                        a <- promoters(g, upstream = abs(min(.ele)), 
                                       downstream = 0)
                        b <- promoters(g, upstream = abs(max(.ele)),
                                       downstream = 0)
                        filterByOverlaps(a, b,
                                         ignore.strand = ignore.strand)
                      },
                      "ups"={
                        promoters(g, upstream = abs(min(.ele)), 
                                  downstream = abs(max(.ele)))
                      },
                      "pp"={
                        a <- promoters(g, upstream = 0, 
                                       downstream = max(.ele))
                        b <- promoters(g, upstream = 0,
                                       downstream = min(.ele))
                        filterByOverlaps(a, b,
                                         ignore.strand = ignore.strand)
                      },
                      "dd"={
                        a <- downstreams(g, upstream = 0, 
                                       downstream = max(.ele))
                        b <- downstreams(g, upstream = 0,
                                       downstream = min(.ele))
                        filterByOverlaps(a, b,
                                         ignore.strand = ignore.strand)
                      },
                      "dws"={
                        downstreams(g, upstream = abs(min(.ele)),
                                    downstream = abs(max(.ele)))
                      },
                      "ww"={
                        a <- downstreams(g, upstream = abs(min(.ele)), 
                                         downstream = 0)
                        b <- downstreams(g, upstream = abs(max(.ele)),
                                         downstream = 0)
                        filterByOverlaps(a, b,
                                         ignore.strand = ignore.strand)
                      })
        if(fil==1){
          .ele <- filterByOverlaps(.ele, g, ignore.strand = ignore.strand)
        }
        seqlevelsStyle(.ele) <- seql[1]
        return(.ele)
      }
      if(is.function(.ele)){
        .ele <- .ele(TxDb)
        if(is(.ele, 'GRangesList')) .ele <- unlist(.ele)
        seqlevelsStyle(.ele) <- seql[1]
        return(.ele)
      }
    })
    l <- lengths(anno)
    n <- names(anno)
    anno <- unlist(GRangesList(anno))
    mcols(anno) <- DataFrame(type=rep(n, l))
    
    ## filter peaks by seqlev
    if(!missing(seqlev)){
      if(length(seqlev)>0){
        peaks <- lapply(peaks, 
                        function(.ele) .ele[seqnames(.ele) %in% seqlev])
      }
    }
    
    peaks <- lapply(peaks, FUN = function(.peaks){
      .peaks$id <- seq_along(.peaks)
      ol <- findOverlaps(.peaks, anno, ignore.strand=ignore.strand)
      .peaks$annoType <- "undefined"
      .p1 <- .peaks[-unique(queryHits(ol))]
      .p2 <- .peaks[queryHits(ol)]
      .p2$annoType <- anno$type[subjectHits(ol)]
      .peaks <- c(.p1, .p2)
      .peaks[order(.peaks$id)]
    })
    if(isGRanges){
      peaks <- peaks[[1]]
    }
    melt <- function(.ele){
      .i <- unique(.ele$id)
      .j <- c(names(breaks), "undefined")
      x <- matrix(0, nrow = length(.i), ncol = length(.j))
      colnames(x) <- .j
      rownames(x) <- .i
      .ele <- split(.ele$id, .ele$annoType)
      for(i in seq_along(.ele)){
        x[.ele[[i]], names(.ele)[i]] <- 1
      }
      as.data.frame(x)
    }
    if(isGRanges){## for upset
      dat <- melt(peaks)
    }else{## bar-plot
      dat <- lapply(peaks, melt)
    }
    
    return(list(peaks=peaks, plotData=dat))
  }
