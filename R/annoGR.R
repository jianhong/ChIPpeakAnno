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
                 re <- "the names has duplicates"
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
              gr <- if(!missing(OrganismDb)) TxDb2GR(ranges, feature, OrganismDb) else TxDb2GR(ranges, feature)
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
              gr <- EnsDb2GR(ranges, feature)
              newAGR(gr=gr, source=source,
                     date=date, feature=feature, 
                     mdata=mdata)
          })
