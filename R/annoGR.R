## Annotation class
## need DBI
setClass("annoGR", 
         representation(source="character",
                        date="Date",
                        feature="character",
                        mdata="data.frame"),
         contains="GRanges",
         validity=function(object){
             re <- TRUE
             if(length(object@seqnames)<1) re <- "annotation is empty"
             if(is.null(names(object@ranges))) 
                 re <- "annotation must have names"
             if(any(duplicated(names(object@ranges))))
                 re <- "the names has dupoicates"
             if(!is.null(object@mdata)){
                 if(!all(colnames(object@mdata)==c("name", "value"))){
                     re <- "colnames of mdata must be name and value"
                 }
             }
             re
         })

newAnnoGR <- function (seqnames = Rle(), 
                       ranges = IRanges(), 
                       strand = Rle("*", length(seqnames)), 
                       mcols = DataFrame(), 
                       seqlengths = NULL, 
                       seqinfo = NULL,
                       ...) 
{
    if (!is(seqnames, "Rle")) 
        seqnames <- Rle(seqnames)
    if (!is.factor(runValue(seqnames))) 
        runValue(seqnames) <- factor(runValue(seqnames), 
                                     levels = unique(runValue(seqnames)))
    if (class(ranges) != "IRanges") 
        ranges <- as(ranges, "IRanges")
    if (!is(strand, "Rle")) 
        strand <- Rle(strand)
    if (!is.factor(runValue(strand)) || !identical(levels(runValue(strand)), 
                                                   levels(strand()))) 
        runValue(strand) <- strand(runValue(strand))
    if (anyMissing(runValue(strand))) {
        warning("missing values in strand converted to \"*\"")
        runValue(strand)[is.na(runValue(strand))] <- "*"
    }
    lx <- max(length(seqnames), length(ranges), length(strand))
    if (lx > 1) {
        if (length(seqnames) == 1) 
            seqnames <- rep(seqnames, lx)
        if (length(ranges) == 1) 
            ranges <- rep(ranges, lx)
        if (length(strand) == 1) 
            strand <- rep(strand, lx)
    }
    if (is.null(seqlengths)) 
        seqlengths <- setNames(rep(NA_integer_, length(levels(seqnames))), 
                               levels(seqnames)) ## stats
    if (is.null(seqinfo)) 
        seqinfo <- Seqinfo(names(seqlengths), seqlengths) ##GenomicInfoDb
    runValue(seqnames) <- factor(runValue(seqnames), seqnames(seqinfo))
    if (!is(mcols, "DataFrame")) 
        stop("'mcols' must be a DataFrame object")
    if (ncol(mcols) == 0L) {
        mcols <- mcols(ranges)
        if (is.null(mcols)) 
            mcols <- new("DataFrame", nrows = length(seqnames))
    }
    if (!is.null(mcols(ranges))) 
        mcols(ranges) <- NULL
    if (!is.null(rownames(mcols))) {
        if (is.null(names(ranges))) 
            names(ranges) <- rownames(mcols)
        rownames(mcols) <- NULL
    }
    new("annoGR", seqnames = seqnames, ranges = ranges, strand = strand, 
        seqinfo = seqinfo, elementMetadata = mcols, ...)
}


newAGR <- function(gr, ...){
    newAnnoGR(seqnames = seqnames(gr), 
              ranges = ranges(gr), 
              strand = strand(gr), 
              mcols = mcols(gr), 
              seqlengths = seqlengths(gr), 
              seqinfo = seqinfo(gr),
              ...)
}

if(!isGeneric("annoGR")){
    setGeneric("annoGR", function(ranges, ...) standardGeneric("annoGR"))
}
if(!isGeneric("info")){
    setGeneric("info", function(object) standardGeneric("info"))
}

setAs(from="annoGR", to="GRanges", function(from){
    do.call(GRanges, args=append(list(seqnames=seqnames(from), 
                                      ranges=ranges(from),
                                      strand=strand(from),
                                      seqlengths=seqlengths(from),
                                      seqinfo=seqinfo(from)),
                                 as.list(mcols(from))))
})
setAs(from="GRanges", to="annoGR", function(from){
    annoGR(from)
})

setMethod("info", "annoGR", function(object){
    cat(class(object), "object;\n")
    cat("# source: ", object@source, "\n")
    cat("# create at: ", format(object@date, "%a %b %d %X %Y %Z"), "\n")
    cat("# feature: ", object@feature, "\n")
    mdata <- object@mdata
    for (i in seq_len(nrow(mdata))) {
        cat("# ", mdata[i, "name"], ": ", mdata[i, "value"],
            "\n", sep="")
    }
})

setMethod("annoGR", "GRanges", 
          function(ranges, feature="group", date, ...){
              if(missing("date")) date <- Sys.Date()
              if(is.null(names(ranges))){
                  names(ranges) <- make.names(
                      formatC(1:length(ranges),
                              width=nchar(length(ranges)),
                              flag="0"))
              }
              newAGR(gr=ranges,
                     date=date, feature=feature,
                     ...)
          })

setMethod("annoGR", "TxDb", 
          function(ranges, feature=c("gene", "transcript", "exon",
                                     "CDS", "fiveUTR", "threeUTR",
                                     "microRNA", "tRNAs", "geneModel"),
                   date, source, mdata, OrganismDb){
              feature <- match.arg(feature)
              if(missing(mdata)) {
                  mdata <- 
                      dbGetQuery(dbconn(ranges), "select * from metadata")
              }
              if(missing(source)) 
                  source <- deparse(substitute(ranges, env=parent.frame()))
              if(missing(date)) date <- Sys.Date()
              gr <- 
                  switch(feature,
                         geneModel={
                             exon <- exonsBy(ranges, "tx", use.names=TRUE)
                             tids <- rep(names(exon), elementNROWS(exon))
                             exon <- unlist(exon)
                             if(length(exon)){
                                 exon$tx_name <- tids
                                 exon$feature_type <- "ncRNA"
                                 cds <- cdsBy(ranges, "tx", use.names=TRUE)
                                 tids <- 
                                     rep(names(cds), elementNROWS(cds))
                                 cds <- unlist(cds)
                                 if(length(cds)){
                                     mcols(cds) <- NULL
                                     cds$tx_name <- tids
                                     cds$feature_type <- "CDS"
                                 }
                                 utr5 <- 
                                     fiveUTRsByTranscript(ranges,
                                                          use.names=TRUE)
                                 tids <- rep(names(utr5), 
                                             elementNROWS(utr5))
                                 utr5 <- unlist(utr5)
                                 if(length(utr5)){
                                     mcols(utr5) <- NULL
                                     utr5$tx_name <- tids
                                     utr5$feature_type <- "5UTR"
                                 }
                                 utr3 <- 
                                     threeUTRsByTranscript(ranges,
                                                           use.names=TRUE)
                                 tids <- rep(names(utr3), 
                                             elementNROWS(utr3))
                                 utr3 <- unlist(utr3)
                                 if(length(utr3)){
                                     mcols(utr3) <- NULL
                                     utr3$tx_name <- tids
                                     utr3$feature_type <- "3UTR"
                                 }
                                 anno <- c(cds, utr5, utr3)
                                 left <- exon[!(exon$tx_name %in% anno$tx_name)]
                                 ##check logical, anno covered all annotation
                                 right <- exon[exon$tx_name %in% anno$tx_name]
                                 rd1 <- reduce(anno)
                                 rd2 <- reduce(right)
                                 if(!identical(rd1, rd2)){
                                     stop("some annotation is missing! bug!")
                                 }
                                 mcols(exon) <- 
                                     mcols(exon)[, 
                                                 c("tx_name", "feature_type")]
                                 exon <- c(exon, anno) ## merge ncRNA with anno
                                 if(!missing(OrganismDb)){
                                     if(class(OrganismDb)=="OrganismDb"){
                                         symbol <- tryCatch(
                                             select(OrganismDb, 
                                                    keys=unique(exon$tx_name),
                                                    columns="SYMBOL",
                                                    keytype="TXNAME"),
                                             error=NULL)
                                         if(length(symbol)>0){
                                             exon$symbol <- 
                                                 symbol[match(exon$tx_name, 
                                                              symbol[, 1]),
                                                        "SYMBOL"]
                                         }
                                     }else{
                                         message("OrganismDb must be an object
                                                 of OrganismDb.")
                                     }
                                     }
                                 ## sort exon
                                 exon <- exon[order(exon$tx_name)]
                                 ### get each tx_name first start pos
                                 tids <- rle(exon$tx_name)
                                 tids$values <- 
                                     tapply(start(exon), exon$tx_name, min)
                                 tids <- inverse.rle(tids)
                                 exon <- 
                                     exon[order(as.character(seqnames(exon)),
                                                tids, 
                                                start(exon))]
                                 names(exon) <- make.names(names(exon), 
                                                           unique=TRUE)
                                 }
                             exon
                         },
                         gene={
                             g <- genes(ranges, columns="gene_id")
                             names(g) <- g$gene_id
                             g$gene_id <- NULL
                             g
                         },
                         exon={
                             e <- exons(ranges, 
                                        columns=c("exon_id", 
                                                  "tx_name", 
                                                  "gene_id"))
                             if(length(e)){
                                 names(e) <- e$exon_id
                                 e$exon_id <- NULL
                             }
                             e
                         },
                         transcript={
                             t <- transcripts(ranges,
                                              columns=c("tx_id",
                                                        "tx_name",
                                                        "gene_id"))
                             if(length(t)){
                                 names(t) <- t$tx_id
                                 t$tx_id <- NULL
                             }
                             t
                         },
                         CDS={
                             c <- cds(ranges, 
                                      columns=c("cds_id",
                                                "tx_name",
                                                "gene_id"))
                             if(length(c)){
                                 names(c) <- c$cds_id
                                 c$cds_id <- NULL
                             }
                             c
                         },
                         fiveUTR={
                             u <- fiveUTRsByTranscript(ranges,
                                                       use.name=TRUE)
                             tids <- rep(names(u), elementNROWS(u))
                             u <- unlist(u)
                             if(length(u)){
                                 u$tx_name <- tids
                                 names(u) <- make.names(names(u), 
                                                        unique=TRUE)
                             }
                             u
                         },
                         threeUTR={
                             u <- threeUTRsByTranscript(ranges,
                                                        use.name=TRUE)
                             tids <- rep(names(u), elementNROWS(u))
                             u <- unlist(u)
                             if(length(u)){
                                 u$tx_name <- tids
                                 names(u) <- make.names(names(u), 
                                                        unique=TRUE)
                             }
                             u
                         },
                         microRNA={
                             m <- microRNAs(ranges)
                             if(length(m)){
                                 names(m) <- m$mirna_id
                                 m$mirna_id <- NULL
                             }
                             m
                         },
                         tRNA=tRNAs(ranges)
                  )
              newAGR(gr=gr, source=source,
                     date=date, feature=feature, 
                     mdata=mdata)
          })

setMethod("annoGR", "EnsDb",
          function(ranges, 
                   feature=c("gene", "transcript", "exon", "disjointExons"),
                   date, source, mdata){
              feature <- match.arg(feature)
              if(missing(mdata)) {
                  mdata <- 
                      dbGetQuery(dbconn(ranges), "select * from metadata")
              }
              if(missing(source)) 
                  source <- deparse(substitute(ranges, env=parent.frame()))
              if(missing(date)) date <- Sys.Date()
              gr <- 
                  switch(feature,
                         disjointExons={
                             e <- disjointExons(ranges,
                                                aggregateGenes=FALSE,
                                                includeTranscripts=TRUE)
                             l <- length(e)
                             names(e) <- make.names(
                                 formatC(seq_len(l), 
                                         width=nchar(as.character(l)),
                                         flag="0"))
                             e
                         },
                         gene={
                             g <- genes(ranges, columns=c("gene_id",
                                                          "gene_name"))
                             names(g) <- g$gene_id
                             g$gene_id <- NULL
                             g
                         },
                         exon={
                             e <- exons(ranges, 
                                        columns=c("exon_id", 
                                                  "tx_id", 
                                                  "gene_id",
                                                  "gene_name"))
                             if(length(e)){
                                 names(e) <- make.names(names(e), 
                                                        unique=TRUE,
                                                        allow_=TRUE)
                             }
                             e
                         },
                         transcript={
                             t <- transcripts(ranges,
                                              columns=c("tx_id",
                                                        "gene_id",
                                                        "gene_name"))
                             if(length(t)){
                                 names(t) <- make.names(names(t), 
                                                        unique=TRUE,
                                                        allow_=TRUE)
                             }
                             t
                         }
                  )
              seqlevelsStyle(gr) <- "UCSC"
              newAGR(gr=gr, source=source,
                     date=date, feature=feature, 
                     mdata=mdata)
          })
