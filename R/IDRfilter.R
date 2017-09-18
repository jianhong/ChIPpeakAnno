IDRfilter <- function(peaksA, peaksB, bamfileA, bamfileB, 
                      maxgap=-1L, minoverlap=0L, singleEnd=TRUE,
                      IDRcutoff=0.01, ...){
    stopifnot(class(peaksA)=="GRanges")
    stopifnot(class(peaksB)=="GRanges")
    stopifnot(file.exists(bamfileA))
    stopifnot(file.exists(bamfileB))
    ol <- findOverlapsOfPeaks(peaksA, peaksB, 
                              maxgap=maxgap, 
                              minoverlap=minoverlap)
    ol <- ol$peaklist[grepl("\\/\\/\\/", names(ol$peaklist))][[1]]
    if(length(ol)<1) return(GRanges())
    names(ol) <- paste0("olp", 
                        formatC(1:length(ol), 
                                width=nchar(as.character(length(ol))),
                                flag="0"))
    coverage <- summarizeOverlaps(features = ol, 
                                  reads = c(bamfileA, bamfileB),
                                  mode=Union, 
                                  ignore.strand = FALSE, 
                                  singleEnd=singleEnd)
    idr <- est.IDR(assay(coverage)/width(rowRanges(coverage)),
                   mu=2.07, sigma=1.34, rho=0.89, p=0.84)
    ol[idr$IDR<IDRcutoff]
}
