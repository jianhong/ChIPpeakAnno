#' @importFrom stats rexp
randPeaks <- function(A, grs, N, ...){
    len <- sapply(grs, length)
    if(any(len<1)){
        stop("All length of grs must be greater than 0")
    }
    if(length(N)!=length(grs)){
        stop("length of N and grs are not identical")
    }
    
    s <- mapply(runifGR, grs, N)
    s <- unlist(GRangesList(s))
    
    wid <- width(A)
    ## how to avoid warning?
    wid <- rexp(length(s),1/(mean(wid)-min(wid))) + min(wid)
    halfwid <- floor(wid/2)
    suppressWarnings(start(s) <- start(s)-halfwid)
    suppressWarnings(width(s) <- 2*halfwid)
    s <- trim(s)
    s
}

runifGR <- function(grNoN, n){
    grNoN <- grNoN[sample.int(length(grNoN))]
    wid <- width(grNoN)
    totL <- sum(as.numeric(wid))
    ends <- cumsum(as.numeric(wid))
    starts <- c(1, ends[-length(ends)] + 1)
    idx <- sample.int(n=totL, size=n, replace=TRUE)
    pos <- cut(idx, breaks=c(0, ends), labels=1:length(grNoN))
    pos <- as.integer(as.character(pos))
    off <- idx - starts[pos] - 1
    gr <- grNoN[pos]
    suppressWarnings(gr <- shift(gr, shift=off))
    gr <- trim(gr)
    width(gr) <- 1
    gr
}