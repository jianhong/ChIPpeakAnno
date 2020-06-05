#' Add metadata of the GRanges objects used for findOverlapsOfPeaks
#' 
#' Add metadata to to overlapping peaks after calling findOverlapsOfPeaks.
#' 
#' 
#' @param ol An object of overlappingPeaks, which is output of
#' \link{findOverlapsOfPeaks}.
#' @param colNames Names of metadata column to be added. If it is NULL,
#' addMetadata will guess what to add.
#' @param FUN A function to be called
#' @param ...  Arguments to the function call.
#' @return return value is An object of \link{overlappingPeaks}.
#' @export
#' @importFrom S4Vectors mcols aggregate
#' @author Jianhong Ou
#' @seealso See Also as \code{\link{findOverlapsOfPeaks}}
#' @keywords misc
#' @examples
#' 
#' peaks1 <- GRanges(seqnames=c(6,6,6,6,5),
#'                  IRanges(start=c(1543200,1557200,1563000,1569800,167889600),
#'                          end=c(1555199,1560599,1565199,1573799,167893599),
#'                          names=c("p1","p2","p3","p4","p5")),
#'                  strand="+",
#'                  score=1:5, id=letters[1:5])
#' peaks2 <- GRanges(seqnames=c(6,6,6,6,5),
#'                   IRanges(start=c(1549800,1554400,1565000,1569400,167888600),
#'                           end=c(1550599,1560799,1565399,1571199,167888999),
#'                           names=c("f1","f2","f3","f4","f5")),
#'                   strand="+",
#'                   score=6:10, id=LETTERS[1:5])
#' ol <- findOverlapsOfPeaks(peaks1, peaks2)
#' addMetadata(ol)
#' 
addMetadata <- function(ol, colNames=NULL, FUN=c, ...){
    stopifnot(class(ol)=="overlappingPeaks")
    stopifnot(length(ol$all.peaks)>0)
    PeaksList <- ol$all.peaks
    if(length(colNames)==0){
        colNames <- lapply(PeaksList, function(.ele) colnames(mcols(.ele)))
        colNames <- table(unlist(colNames))
        colNames <- names(colNames)[colNames==length(PeaksList)]
    }
    stopifnot(is.character(colNames))
    stopifnot(length(colNames)>0)
    stopifnot(is.function(FUN))
    for(i in 1:length(PeaksList)){
        if(!all(colNames %in% colnames(mcols(PeaksList[[i]])))){
            stop(paste("colNames: ", 
                       paste(colNames[!colNames %in% 
                                          colnames(mcols(PeaksList[[i]]))], 
                             collapse=", "),
                       "does not exist in the metadata of all the list."))
        }
    }
    PeaksList <- lapply(PeaksList, function(.ele) .ele[, colNames])
    cl <- lapply(PeaksList, function(.ele) sapply(mcols(.ele), class))
    for(i in 2:length(cl)){
        if(!identical(cl[[1]], cl[[i]])){
            stop("class of metadata are not identical.")
        }
    }
    Peaks <- unlist(GRangesList(PeaksList), use.names = FALSE)
    ol$peaklist <- lapply(ol$peaklist, function(.ele){
        .mcols <- aggregate(mcols(Peaks[unlist(.ele$peakNames)]), 
                            by=list(
                                addMetatdata_group=
                                    rep(seq_along(.ele), 
                                        lengths(.ele$peakNames))), 
                            FUN=mean, ...)
        .mcols <- .mcols[order(.mcols$addMetatdata_group), ]
        .n <- ncol(mcols(.ele))
        mcols(.ele) <- cbind(mcols(.ele), .mcols[, colNames])
        if(.n>0){
          colnames(mcols(.ele))[-seq.int(.n)] <- colNames
        }else{
          colnames(mcols(.ele)) <- colNames
        }
        .ele
    })
    ol
}
