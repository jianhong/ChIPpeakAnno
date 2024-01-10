#' plot enrichment results
#' 
#' Plot the GO/KEGG/reactome enrichment results
#' 
#' @param res output of \link{getEnrichedGO}, \link{getEnrichedPATH}.
#' @param n number of terms to be plot.
#' @param style plot vertically or horizontally
#' @param strlength shorten the description of term by the number of char.
#' @param orderBy order the data by pvalue, termId or none.
#' @param label_wrap soft wrap the labels (i.e. descriptions of the GO or PATHWAY terms), default to 40 characters.
#' @param label_substring_to_remove remove common substring from label, 
#' default to NULL. Special characters must be escaped. E.g. if you would like 
#' to remove "Homo sapiens (human)" from labels, you must use "Homo sapiens \\\(
#' human\\\)".
#' @author Jianhong Ou, Kai Hu
#' @return an object of ggplot
#' @importFrom ggplot2 ggplot aes geom_bar geom_point scale_x_discrete scale_y_continuous geom_text
#' xlab ylab theme_classic theme facet_grid expansion element_text
#' @importFrom stats as.formula
#' @importFrom scales label_wrap
#' @export
#' @examples 
#' data(enrichedGO)
#' enrichmentPlot(enrichedGO)
#' if (interactive()||Sys.getenv("USER")=="jianhongou") {
#'      
#'      library(org.Hs.eg.db)
#'      library(GO.db)
#'      bed <- system.file("extdata", "MACS_output.bed", package="ChIPpeakAnno")
#'      gr1 <- toGRanges(bed, format="BED", header=FALSE) 
#'      gff <- system.file("extdata", "GFF_peaks.gff", package="ChIPpeakAnno")
#'      gr2 <- toGRanges(gff, format="GFF", header=FALSE, skip=3)
#'      library(EnsDb.Hsapiens.v75) ##(hg19)
#'      annoData <- toGRanges(EnsDb.Hsapiens.v75)
#'      gr1.anno <- annoPeaks(gr1, annoData)
#'      gr2.anno <- annoPeaks(gr2, annoData)
#'      over <- lapply(GRangesList(gr1=gr1.anno, gr2=gr2.anno), 
#'                     getEnrichedGO, orgAnn="org.Hs.eg.db",
#'                     maxP=.05, minGOterm=10, condense=TRUE)
#'      enrichmentPlot(over$gr1)
#'      enrichmentPlot(over$gr2, style = "h")
#'  }
enrichmentPlot <- function(res, n=20, strlength=Inf,
                           style = c("v", "h"),
                           label_wrap = 40,
                           label_substring_to_remove = NULL,
                           orderBy=c("pvalue", "termId", "none")){
  if(is.data.frame(res)){## output of getEnrichedPATH
    res <- list(path=res)
  }
  stopifnot("n must be a numeric(1)"=is.numeric(n) && length(n)==1)
  orderBy <- match.arg(orderBy)
  style <- match.arg(style)
  stopifnot(is.integer(as.integer(label_wrap)))
  if(is.list(res[[1]])&&is.data.frame(res[[1]][[1]])){
    ## list of list, output of getEnrichedGO for multiple samples
    ## dot plot
    res <- swapList(res)
    res <- lapply(res, function(.ele){
      ## unlist and add source
      .e <- do.call(rbind, .ele)
      .e$source <- rep(names(.ele), vapply(.ele, nrow, FUN.VALUE = 0))
      .e
    })
  }
  p <- lapply(res, function(.ele){
    cn <- colnames(.ele)
    cn.id <- cn[grepl("\\.id$", cn)]
    cn.term <- cn[grepl("\\.term$", cn)]
    if(length(.ele$source)!=nrow(.ele)){
      .ele$source <- "undefined"
    }
    if(nrow(.ele)==0) return(data.frame())
    if(!all(c(cn.id, cn.term, "pvalue", "count.InDataset", 
             "count.InGenome", "source") %in% colnames(.ele))){
      return(data.frame())
    }
    plotdata <- .ele[!is.na(.ele$pvalue), 
                     c(cn.id, cn.term, "pvalue", "count.InDataset", 
                       "count.InGenome", "source")]
    plotdata <- as.data.frame(plotdata)
    plotdata <- unique(plotdata)
    plotdata$qvalue <- -1*log10(plotdata$pvalue)
    plotdata <- switch(orderBy,
                       "pvalue"=plotdata[order(plotdata$pvalue), ],
                       "termId"=plotdata[order(plotdata[, cn.term]), ],
                       plotdata)
    if(nrow(plotdata)>n) plotdata <- plotdata[seq.int(n), ]
    plotdata$Description <- shortStrs(plotdata[, cn.term], len=strlength)
    plotdata$Count <- plotdata$count.InDataset
    plotdata$GeneRatio <- plotdata$count.InDataset/plotdata$count.InGenome
    plotdata
  })
  plotdata <- do.call(rbind, p)
  plotdata$category <- rep(names(res), vapply(p, nrow, FUN.VALUE = 0))
  plotdata$Description <- removeLabelSubstring(plotdata, label_substring_to_remove)

  if(nrow(plotdata)<2){
    warning("two less data to plot")
    return(plotdata)
  }

  if(all(plotdata$source=="undefined")){
    if (style == "v") {
      p <- ggplot(plotdata, 
                  aes(x = switch(orderBy,
                                 "pvalue" = reorder(Description, pvalue),
                                 "termId" = reorder(Description, plotdata[, 2]),
                                 "Description"),
                      y = qvalue, fill = Count, label = Count))
    } else if (style == "h") {
      p <- ggplot(plotdata, 
                  aes(x = switch(orderBy,
                                 "pvalue" = reorder(Description, pvalue, decreasing = TRUE),
                                 "termId" = reorder(Description, plotdata[, 2]),
                                 "Description"),
                      y = qvalue, fill = Count, label = Count))
    }
    p <- p +
      geom_bar(stat="identity") +
      scale_x_discrete(label = label_wrap(as.integer(label_wrap))) +
      scale_y_continuous(expand = expansion(mult = c(0, .1))) +
      xlab("") + ylab("-log10(p-value)") +
      theme_classic() + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + 
      facet_grid(~ category, scales = "free_x", space = "free_x")
    
    if (style == "v") {
      p <- p + 
            geom_text(vjust=-.1)
    } else if (style == "h") {
      p <- p + 
            geom_text(hjust=-.1) + 
            coord_flip()
    }
    p
  }else{
    ## multiple samples dot plot
    if (style == "v") {
      p <- ggplot(plotdata,
                  aes(y = switch(orderBy, 
                                 "pvalue" = reorder(Description, pvalue),
                                 "termId" = reorder(Description, plotdata[, 2]),
                                 Description),
                      x = source, color = qvalue, size = GeneRatio))
    } else if (style == "h") {
      p <- ggplot(plotdata, 
                  aes(y = switch(orderBy,
                                 "pvalue" = reorder(Description, pvalue, decreasing = TRUE),
                                 "termId" = reorder(Description, plotdata[, 2]),
                                 Description),
                      x = source, color = qvalue, size = GeneRatio))
    }
    p <- p +
      geom_point() + theme_classic() +
      facet_grid(~ category, scales = "free_x", space = "free_x")
  
    if (style == "v") {
      p
    } else if (style == "h") {
      p <- p + coord_flip()
      p
    }
  }
}

removeLabelSubstring <- function(plotdata, label_substring_to_remove = NULL) {
  if (!is.null(label_substring_to_remove)) {
    gsub(label_substring_to_remove, "", plotdata$Description)
  } else {
    plotdata$Description
  }
}

shortStrs <- function(strs, len=60){
  if(length(strs)==0) return(strs)
  strs <- as.character(strs)
  shortStr <- function(str, len=60){
    stopifnot(length(str)==1)
    stopifnot(is.character(str))
    if(nchar(str)<=len) return(str)
    strs <- strsplit(str, " ")[[1]]
    nc <- nchar(strs)
    nclast <- nc[length(nc)] + 3
    paste0(substring(str, first = 1, last = len-nclast), "...",
           strs[length(strs)])
  }
  strs <- sapply(strs, shortStr, len=len)
  make.unique(strs)
}
