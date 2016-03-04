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
        .mcols <- sapply(colNames, function(i){
            sapply(.ele$peakNames, function(.name) FUN(mcols(Peaks[.name])[, i], ...))
        }, simplify = FALSE)
        mcols(.ele) <- cbind(mcols(.ele), DataFrame(as.data.frame(do.call(cbind, .mcols), stringsAsFactors=FALSE)))
        .ele
    })
    ol
}