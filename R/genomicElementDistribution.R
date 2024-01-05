#' Genomic Element distribution
#' 
#' Plot pie chart for genomic element distribution
#' 
#' @details The distribution will be calculated by geneLevel,
#' ExonIntron, and Exons The geneLevel will be categorized as
#' promoter region, gene body, gene downstream and distal intergenic region.
#' The ExonIntron will be categorized as exon, intron and intergenic.
#' The Exons will be categorized as 5' UTR, 3'UTR and CDS.
#' The precedence will follow the order of labels defination.
#' For example, for ExonIntron, if a peak overlap with both exon and intron, 
#' and exon is specified before intron, then only exon will
#' be incremented for the same example.
#' @param peaks peak list, \link[GenomicRanges:GRanges-class]{GRanges} object or
#' a \link[GenomicRanges:GRangesList-class]{GRangesList}.
#' @param TxDb an object of \code{\link[GenomicFeatures:TxDb-class]{TxDb}}
#' @param seqlev sequence level should be involved. 
#' Default is all the sequence levels in intersect of peaks and TxDb.
#' @param nucleotideLevel Logical. Choose between peak centric and nucleotide
#' centric view. Default=FALSE
#' @param ignore.strand logical. Whether the strand of the input ranges
#' should be ignored or not. Default=TRUE
#' @param promoterRegion numeric. The upstream and downstream of genes to define
#' promoter region.
#' @param promoterLevel list. The breaks, labels, and colors for divided range 
#' of promoters. The breaks must be from 5' -> 3' and the percentage will use
#' the fixed precedence 3' -> 5' 
#' @param geneDownstream numeric. The upstream and downstream of genes to define
#' gene downstream region.
#' @param labels list. A list for labels for the genomic elements.
#' @param labelColors named character vector. The colors for each labels.
#' @param plot logic. Plot the pie chart for the genomic elements or not.
#' @param keepExonsInGenesOnly logic. Keep the exons within annotated gene only.
#' @export
#' @importFrom ggplot2 ggplot geom_rect xlim coord_polar aes_string geom_bar
#' coord_flip scale_fill_manual theme_void theme_bw facet_wrap geom_col 
#' geom_text guide_legend
#' @importFrom GenomicFeatures intronsByTranscript exons fiveUTRsByTranscript 
#' threeUTRsByTranscript genes
#' @importFrom stats as.formula
#' @return Invisible list of data for plot.
#' @examples 
#' if (interactive() || Sys.getenv("USER")=="jianhongou"){
#'   data(myPeakList)
#'   if(require(TxDb.Hsapiens.UCSC.hg19.knownGene)){
#'   seqinfo(myPeakList) <- 
#'   seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene)[seqlevels(myPeakList)]
#'   myPeakList <- GenomicRanges::trim(myPeakList)
#'   myPeakList <- myPeakList[width(myPeakList)>0]
#'     genomicElementDistribution(myPeakList, 
#'         TxDb.Hsapiens.UCSC.hg19.knownGene)
#'     genomicElementDistribution(myPeakList, 
#'         TxDb.Hsapiens.UCSC.hg19.knownGene,
#'         nucleotideLevel = TRUE)
#'     genomicElementDistribution(myPeakList, 
#'         TxDb.Hsapiens.UCSC.hg19.knownGene,
#'         promoterLevel=list(
#'         #from 5' -> 3', fixed precedence 3' -> 5'
#'         breaks = c(-2000, -1000, -500, 0, 100),
#'         labels = c("upstream 1-2Kb", "upstream 0.5-1Kb", 
#'                    "upstream <500b", "TSS - 100b"),
#'         colors = c("#FFE5CC", "#FFCA99", 
#'                    "#FFAD65", "#FF8E32")))
#'   }
#' }
genomicElementDistribution <- 
  function(peaks, TxDb, seqlev, nucleotideLevel=FALSE, ignore.strand=TRUE,
           promoterRegion=c(upstream=2000, downstream=100),
           geneDownstream=c(upstream=0, downstream=1000),
           labels=list(geneLevel=c(promoter="Promoter",
                                   geneDownstream="Downstream",
                                   geneBody="Gene body",
                                   distalIntergenic="Distal Intergenic"),
                       ExonIntron=c(exon="Exon",
                                    intron="Intron",
                                    intergenic="Intergenic"),
                       Exons=c(utr5="5' UTR",
                               utr3="3' UTR",
                               CDS="CDS",
                               otherExon="Other exon"),
                       group=c(geneLevel="Transcript Level",
                               promoterLevel="Promoter Level",
                               Exons="Exon level",
                               ExonIntron="Exon/Intron/Intergenic")),
           labelColors = c(promoter="#D3FF00",
                           geneBody="#9EFF00",
                           geneDownstream="#6AFF00",
                           distalIntergenic="#00FF28",
                           exon="#6600FF",
                           intron="#8F00FF",
                           intergenic="#DA00FF",
                           utr5="#00FFDB",
                           utr3="#00DFFF",
                           CDS="#00A0FF",
                           otherExon="#006FFF"),
           plot = TRUE,
           keepExonsInGenesOnly = TRUE,
           promoterLevel){
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
    stopifnot("nuleotideLevel is not logical"=is.logical(nucleotideLevel))
    stopifnot("promoterRegion should contain element upstream and downstream"=
                all(c("upstream", "downstream") %in% names(promoterRegion)))
    stopifnot("geneDownstream should contain element upstream and downstream"=
                all(c("upstream", "downstream") %in% names(geneDownstream)))
    stopifnot("Elements in promoterRegion should be numeric"=
                is.numeric(promoterRegion))
    stopifnot("Elements in geneDownstream should be numeric"=
                is.numeric(geneDownstream))
    labs <- list(geneLevel=c(promoter="Promoter",
                             geneDownstream="Downstream",
                             geneBody="Gene body",
                             distalIntergenic="Distal Intergenic"),
                 ExonIntron=c(exon="Exon",
                              intron="Intron",
                              intergenic="Intergenic"),
                 Exons=c(utr5="5' UTR",
                         utr3="3' UTR",
                         CDS="CDS",
                         otherExon="Other exon"))
    groupLabels <- c(geneLevel="Gene Level",
                     promoterLevel="Promoter Level",
                     Exons="Exon level",
                     ExonIntron="Exon/Intron/Intergenic")
    labelCols = c(promoter="#D55E00",
                  geneDownstream="#E69F00",
                  geneBody="#51C6E6",
                  distalIntergenic="#AAAAAA",
                  exon="#009DDA",
                  intron="#666666",
                  intergenic="#DDDDDD",
                  utr5="#0072B2",
                  utr3="#56B4E9",
                  CDS="#0033BF",
                  otherExon="#009E73",
                  undefined="#FFFFFF")
    labelCols[names(labelColors)] <- labelColors
    for(i in names(labs)){
      if(i %in% names(labels)){
        ## keep the orders in labels
        n <- intersect(names(labs[[i]]), names(labels[[i]]))
        labs[[i]] <- c(labels[[i]][n], labs[[i]][!names(labs[[i]]) %in% n])
      }
    }
    
    for(i in names(groupLabels)){
      if("group" %in% names(labels)){
        groupLabels[names(labels[["group"]])] <- labels[["group"]]
      }
    }
    if(!missing(promoterLevel)){
      # stopifnot("promoterLevel must be within promoterRegion"=
      #             all(abs(promoterLevel$breaks[promoterLevel$breaks<0])<=
      #                   promoterRegion["upstream"]) && 
      #             all(abs(promoterLevel$breaks[promoterLevel$breaks>0])<=
      #                   promoterRegion["downstream"]))
      promoterLevel$breaks <- sort(promoterLevel$breaks)
      stopifnot("breaks, labels and colors of promoterLevel are not paired"=
                  length(promoterLevel$breaks)==
                  length(promoterLevel$labels)+1 &&
                  length(promoterLevel$labels)==
                  length(promoterLevel$colors))
      proK <- paste0("promoter", seq_along(promoterLevel$labels))
      promoterLevel$upstream <- ifelse(promoterLevel$breaks<0,
                                       abs(promoterLevel$breaks),
                                       0)
      promoterLevel$upstream <- 
        promoterLevel$upstream[-length(promoterLevel$upstream)]
      promoterLevel$downstream <- ifelse(promoterLevel$breaks>0,
                                         promoterLevel$breaks,
                                         0)
      promoterLevel$downstream <- promoterLevel$downstream[-1]
      proV <- promoterLevel$labels
      names(proV) <- proK
      labs <- c(list("promoterLevel"=proV), 
                labs)
      proV <- promoterLevel$colors
      names(proV) <- proK
      labelCols <- c(proV, labelCols)
    }else{
      promoterLevel <- NULL
    }
    
    defaultW <- getOption("warn")
    options(warn = -1)
    on.exit({warn = defaultW})
    
    seql <- seqlevelsStyle(peaks)
    
    ## set annotation
    anno <- GRangesList()
    for(i in names(labs)){
      anno[[i]] <- switch (i,
        "promoterLevel" = {
          suppressMessages(g <- genes(TxDb, single.strand.genes.only=TRUE))
          pro <- mapply(FUN=function(upstream, downstream){
            promoters(g, upstream = upstream, downstream = downstream)
          }, promoterLevel$upstream, 
          promoterLevel$downstream,
          SIMPLIFY = FALSE)
          names(pro) <- names(labs[["promoterLevel"]])
          pro <- rev(pro)
          current_anno <- GRanges()
          for(j in seq_along(pro)){
            ca <- filterByOverlaps(pro[[j]], current_anno,
                                   ignore.strand = ignore.strand)
            mcols(ca) <- DataFrame(type=rep(names(pro)[j], length(ca)))
            current_anno <- c(current_anno, ca)
          }
          seqlevelsStyle(current_anno) <- seql[1]
          current_anno
        },
        "geneLevel" = {
          suppressMessages(g <- genes(TxDb, single.strand.genes.only=TRUE))
          pro <- promoters(g, 
                           upstream = promoterRegion["upstream"],
                           downstream = promoterRegion["downstream"])
          dws <- downstreams(g,
                             upstream = geneDownstream["upstream"],
                             downstream = geneDownstream["downstream"])
          pro <- GenomicRanges::trim(pro)
          dws <- GenomicRanges::trim(dws)
          intergenic <- gaps(reduce(c(pro, g, dws), ignore.strand=FALSE))
          intergenic <- intergenic[!strand(intergenic) %in% "*"]
          current_anno <- GRanges()
          ## set precedence
          for(j in names(labs[["geneLevel"]])){
            s <- 
              switch(j,
                     "promoter" =pro,
                     "geneDownstream" =dws,
                     "geneBody" =g,
                     "distalIntergenic" =intergenic)
            ca <- 
              filterByOverlaps(s, current_anno,
                               ignore.strand = ignore.strand)
            mcols(ca) <- DataFrame(type=rep(j, length(ca)))
            current_anno <- c(current_anno, ca)
          }
          seqlevelsStyle(current_anno) <- seql[1]
          current_anno
        },
        "ExonIntron" = {
          exon <- exons(TxDb)
          intron <- unlist(intronsByTranscript(TxDb))
          if(keepExonsInGenesOnly){
            suppressMessages(g <- genes(TxDb, single.strand.genes.only=TRUE))
            ole <- findOverlaps(exon, g, type = "within")
            oli <- findOverlaps(intron, g, type = "within")
            ole <- !seq_along(exon) %in% queryHits(ole)
            oli <- !seq_along(intron) %in% queryHits(oli)
            if(sum(ole)>0 || sum(oli)>0){
              warning(paste(sum(ole), "exons were dropped because there is no",
                            "relative gene level annotations.",
                            sum(oli), "introns were dropped because there is no",
                            "relative gene level annotations."))
              exon <- exon[!ole]
              intron <- intron[!oli]
            }
          }
          intergenic <- gaps(reduce(c(exon, intron), ignore.strand=FALSE))
          intergenic <- intergenic[!strand(intergenic) %in% "*"] 
          current_anno <- GRanges()
          ## set precedence
          for(j in names(labs[["ExonIntron"]])){
            s <- 
              switch(j,
                     "exon" =exon,
                     "intron" =intron,
                     "intergenic" =intergenic)
            ca <- 
              filterByOverlaps(s, current_anno,
                               ignore.strand = ignore.strand)
            mcols(ca) <- DataFrame(type=rep(j, length(ca)))
            current_anno <- c(current_anno, ca)
          }
          seqlevelsStyle(current_anno) <- seql[1]
          current_anno
        },
        "Exons" = {
          utr5 <- unlist(fiveUTRsByTranscript(TxDb))
          utr3 <- unlist(threeUTRsByTranscript(TxDb))
          CDS <- cds(TxDb)
          exon <- exons(TxDb)
          if(keepExonsInGenesOnly){
            suppressMessages(g <- genes(TxDb, single.strand.genes.only=TRUE))
            ole <- findOverlaps(exon, g, type = "within")
            olc <- findOverlaps(CDS, g, type = "within")
            ol5 <- findOverlaps(utr5, g, type = "within")
            ol3 <- findOverlaps(utr3, g, type = "within")
            ole <- !seq_along(exon) %in% queryHits(ole)
            olc <- !seq_along(CDS) %in% queryHits(olc)
            ol5 <- !seq_along(utr5) %in% queryHits(ol5)
            ol3 <- !seq_along(utr3) %in% queryHits(ol3)
            if(sum(ole)>0 || sum(olc)>0 || sum(ol5)>0 || sum(ol3)){
              warning(paste(sum(ole), "exons were dropped because there is no",
                            "relative gene level annotations.",
                            sum(olc), "CDS were dropped because there is no",
                            "relative gene level annotations.",
                            sum(ol5), "utr5 were dropped because there is no",
                            "relative gene level annotations.",
                            sum(ol3), "utr3 were dropped because there is no",
                            "relative gene level annotations."))
              exon <- exon[!ole]
              CDS <- CDS[!olc]
              utr5 <- utr5[!ol5]
              utr3 <- utr3[!ol3]
            }
          }
          current_anno <- GRanges()
          ## set precedence
          for(j in names(labs[["Exons"]])){
            s <- 
              switch(j,
                     "utr5" =utr5,
                     "utr3" =utr3,
                     "CDS" =CDS,
                     "otherExon"=exon)
            ca <- 
              filterByOverlaps(s, current_anno,
                               ignore.strand = ignore.strand)
            mcols(ca) <- DataFrame(type=rep(j, length(ca)))
            current_anno <- c(current_anno, ca)
          }
          seqlevelsStyle(current_anno) <- seql[1]
          current_anno
        }
      )
    }
    groupLabels <- groupLabels[names(anno)]
    groupLabels[is.na(groupLabels)] <- names(anno)[is.na(groupLabels)]
    
    ## filter peaks by seqlev
    if(!missing(seqlev)){
      if(length(seqlev)>0){
        peaks <- lapply(peaks, 
                        function(.ele) .ele[seqnames(.ele) %in% seqlev])
      }
    }
    
    if(nucleotideLevel){
      peaks <- lapply(peaks, function(.ele){
        y <- disjoin(c(.ele, unlist(anno)), ignore.strand=ignore.strand)
        subsetByOverlaps(y, .ele, ignore.strand=ignore.strand)
      })
    }
    
    peaks <- lapply(peaks, FUN = function(.peaks){
      pct <- lapply(anno, FUN = function(.ele){
        y <- .peaks
        ol <- findOverlaps(y, .ele, ignore.strand=ignore.strand)
        ol <- as.data.frame(ol)
        ol <- ol[order(ol$queryHits, ol$subjectHits), ]
        ol <- ol[!duplicated((ol$queryHits)), ]
        y$anno <- rep("undefined", length(y))
        y$anno[ol$queryHits] <- .ele$type[ol$subjectHits]
        y$anno
      })
      
      pct <- do.call(cbind, pct)
      mcols(.peaks) <- cbind(mcols(.peaks), pct)
      .peaks
    })
    if(isGRanges){
      peaks <- peaks[[1]]
    }
    
    melt <- function(m){
      if(nucleotideLevel){
        m <- m[rep(seq_along(m), width(m))]
      }
      m <- mcols(m)[, names(anno)]
      p <- lapply(names(anno), function(.ele){
        tt <- table(m[, .ele])
        data.frame(category=factor(rep(groupLabels[.ele],
                                       length(tt)), 
                                   levels = groupLabels),
                   type=factor(names(tt), levels = rev(names(labelCols))), 
                   percentage=as.numeric(tt)/sum(tt))
      })
      do.call(rbind, p)
    }
    
    if(isGRanges){## donut-plot
      #dat <- reshape2::melt()
      dat <- melt(peaks)
      l <- unlist(unname(labs))
      l1 <- paste0(l, " (", 
                   round(dat$percentage[match(names(l), dat$type)]*100,
                         digits = 1), "%)")
      names(l1) <- names(l)
      l1 <- c(l1, undefined="")
      p <- ggplot(dat, 
                  aes_string(x="category", y="percentage", 
                             fill="type"))+
        geom_col() +
        coord_polar("y") + 
        geom_text(data = subset(dat, !duplicated(dat$category)),
                  aes_string(x="category", label="category"),
                  y=1) +
        theme_void()

    }else{## bar-plot
      dat <- lapply(peaks, melt)
      dat1 <- do.call(rbind, dat)
      dat1$source <- rep(names(peaks), vapply(dat, nrow, FUN.VALUE = 0))
      l1 <- c(unlist(unname(labs)), undefined="")
      p <- ggplot(dat1, 
                  aes_string(x="source", y="percentage", fill="type")) +
        geom_bar(stat="identity") + coord_flip() +
        facet_wrap(as.formula("~ category"), ncol = 1) + 
        theme_bw()
    }
    p <- p + 
      scale_fill_manual(values = labelCols, labels=l1, name=NULL,
                        guide = guide_legend(reverse=TRUE))
    
    if(plot){
      print(p)
    }
    return(invisible(list(peaks=peaks, plot=p)))
  }
