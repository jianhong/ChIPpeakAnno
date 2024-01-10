#' Class \code{annoGR}
#' 
#' An object of class \code{annoGR} represents the annotation data could be
#' used by annotationPeakInBatch.
#' 
#' 
#' @name annoGR-class
#' @rdname annoGR
#' @aliases annoGR-class annoGR
#' @docType class
#' @param ranges an object of \link[GenomicRanges:GRanges-class]{GRanges},
#' \link[GenomicFeatures:TxDb-class]{TxDb} or \link[ensembldb]{EnsDb}
#' @param feature annotation type
#' @param date a \link{Date} object
#' @param ... could be following parameters
#' @param source character, where the annotation comes from
#' @param mdata data frame, metadata from annotation
#' @param OrganismDb an object of OrganismDb. It is used for extracting gene
#' symbol for geneModel group for \link[GenomicFeatures:TxDb-class]{TxDb}
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("annoGR", date, elementMetadata, feature, mdata, ranges, seqinfo,
#' seqnames, source, strand)}
#' @slot seqnames,ranges,strand,elementMetadata,seqinfo slots inherit from 
#' \link[GenomicRanges:GRanges-class]{GRanges}. 
#' The ranges must have unique names.
#' @slot source character, where the annotation comes from
#' @slot date a \link{Date} object
#' @slot feature annotation type, could be "gene", "exon", "transcript", "CDS",
#' "fiveUTR", "threeUTR", "microRNA", "tRNAs", "geneModel" for 
#' \link[GenomicFeatures:TxDb-class]{TxDb} object, or "gene", "exon",
#' "transcript" for \link[ensembldb]{EnsDb} object
#' @slot mdata data frame, metadata from annotation
#' @author Jianhong Ou
#' @keywords classes
#' @exportClass annoGR
#' @import methods
#' @examples
#' 
#'     if(interactive() || Sys.getenv("USER")=="jianhongou"){
#'         library(EnsDb.Hsapiens.v79)
#'         anno <- annoGR(EnsDb.Hsapiens.v79)
#'     }
#' 
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

#' @importFrom S4Vectors Rle DataFrame 
#' @importFrom stats setNames
#' @importFrom GenomeInfoDb Seqinfo
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
    if (!is(ranges, "IRanges")) 
        ranges <- as(ranges, "IRanges")
    if (!is(strand, "Rle")) 
        strand <- Rle(strand)
    if (!is.factor(runValue(strand)) || !identical(levels(runValue(strand)), 
                                                   levels(strand()))) 
        runValue(strand) <- strand(runValue(strand))
    if (any(is.na(runValue(strand)))) {
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
        if (is.null(mcols)){
            rn <- names(ranges)
            removenames <- FALSE
            if(length(rn)!=length(ranges)){
                removenames <- TRUE
                rn <- seq_along(ranges)
            }
            mcols <- DataFrame(row.names = rn)
            if(removenames){
                rownames(mcols) <- NULL
            }
        }
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

setGeneric("annoGR", function(ranges, ...) standardGeneric("annoGR"))
setGeneric("info", function(object) standardGeneric("info"))

#' @name coerce
#' @import GenomicRanges
#' @rdname annoGR
#' @aliases coerce,GRanges,annoGR-method
#' coerce,annoGR,GRanges-method
#' @exportMethod coerce
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

#' @rdname annoGR
#' @exportMethod info
#' @param object annoGR object.
#' @aliases info info,annoGR-method
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

#' @rdname annoGR
#' @exportMethod annoGR
#' @aliases annoGR,GRanges-method
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

#' @importFrom DBI dbGetQuery 
#' @importFrom BiocGenerics dbconn
#' @rdname annoGR
#' @exportMethod annoGR
#' @aliases annoGR,TxDb-method
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

#' @importFrom DBI dbGetQuery 
#' @importFrom BiocGenerics dbconn
#' @rdname annoGR
#' @exportMethod annoGR
#' @aliases annoGR,EnsDb-method
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
