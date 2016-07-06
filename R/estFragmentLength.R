estFragmentLength <- function(bamfiles, index=bamfiles, plot=TRUE, lag.max=1000, ...){
    message("The fragment size is calculated for DNA-seq.")
    res <- mapply(function(f, i) {
        if(suppressMessages(testPairedEndBam(f, index=i))){
            isize <- 
                sort(abs(scanBam(f, index=i, 
                                 param=ScanBamParam(flag=scanBamFlag(isPaired = TRUE, 
                                                                     isProperPair = TRUE,
                                                                     isSecondaryAlignment = FALSE,
                                                                     isFirstMateRead = TRUE,
                                                                     isNotPassingQualityControls = FALSE),
                                                    what="isize"))[[1]]$isize))
            isize <- table(isize)
            isize.x <- as.numeric(names(isize))
            isize.y <- isize
            ## do smooth and select the max point
            isize.smooth <- loess.smooth(isize.x, isize.y, span=1/3, evaluation = length(isize.x))
            pos <- round(isize.smooth$x[which.max(isize.smooth$y)], digits = 0)
            if(plot){
                try({
                    plot(isize.x, isize.y, xlab="Lag", ylab="ACF", main=f, type="l")
                    abline(v=pos, col="red")
                })
            }
            pos
        }else{
            ## count reads for all chromosomes
            seqinfo <- lapply(scanBamHeader(f), 
                                              function(.ele) .ele$targets)
            seqn <- sort(table(unlist(lapply(seqinfo, names))))
            seqn <- names(seqn[which(seqn==max(seqn))])
            seqinfo <- unique(do.call(rbind, lapply(seqinfo, `[`, seqn)))[1, ]
            seqinfo.gr <- GRangesList(lapply(1:length(seqinfo), function(.ele){
                GRanges(names(seqinfo)[.ele], IRanges(1, seqinfo[.ele]))
            }))
            names(seqinfo.gr) <- names(seqinfo)
            cnts <- countBam(f, index=i, 
                             param=ScanBamParam(flag=scanBamFlag(isPaired = FALSE, 
                                                                 isUnmappedQuery = FALSE,
                                                                 isSecondaryAlignment = FALSE,
                                                                 isNotPassingQualityControls = FALSE),
                                                which=seqinfo.gr))
            cnts <- cnts[cnts$nucleotides>0, , drop=FALSE]
            if(nrow(cnts)<1){
                stop("no reads detected.")
            }
            depth <- cnts$nucleotides/cnts$width
            ## select middle one
            select.chr <- as.character(cnts$space[which.min(abs(depth-median(depth)))])
            #read reads
            reads <- readGAlignments(f, index=i, 
                                     param=ScanBamParam(flag=scanBamFlag(isPaired = FALSE, 
                                                                         isUnmappedQuery = FALSE,
                                                                         isSecondaryAlignment = FALSE,
                                                                         isNotPassingQualityControls = FALSE),
                                                        which=seqinfo.gr[[select.chr]],
                                                        what=scanBamWhat()))
            #convert to GRanges
            reads.gr <- as(reads, "GRanges")
            mcols(reads.gr) <- NULL
            reads.gr <- split(reads.gr, strand(reads.gr))
            reads.gr <- sapply(reads.gr, coverage)
            reads.gr.pos <- as.integer(reads.gr[["+"]][[select.chr]])
            reads.gr.neg <- as.integer(reads.gr[["-"]][[select.chr]])
            ## remove 0 from both end
            pos0 <- max(which(reads.gr.pos!=0)[1], which(reads.gr.neg!=0)[1])
            pos1 <- length(reads.gr.pos) - 
                max(which(rev(reads.gr.pos)!=0)[1], which(rev(reads.gr.neg)!=0)[1])
            if(pos1>pos0+5000){
                reads.gr.pos <- reads.gr.pos[pos0:pos1]
                reads.gr.neg <- reads.gr.neg[pos0:pos1]
            }
            if(length(reads.gr.pos)>100000){
                block <- 50000
                len <- ceiling(length(reads.gr.pos)/block)
                signal <- 
                    split(reads.gr.pos, 
                          rep(formatC(1:len, width = nchar(as.character(len)), flag="0"), 
                              each=block)[1:length(reads.gr.pos)])
                signal <- sapply(signal, mean)
                signal <- signal[-len]
                sid <- which.max(signal)
                pos0 <- (sid - 1) * block + 1
                pos1 <- min(sid * block, length(reads.gr.pos))
                reads.gr.pos <- reads.gr.pos[pos0:pos1]
                reads.gr.neg <- reads.gr.neg[pos0:pos1]
            }
            ## cross-correlation 
            reads.gr.pos.ts <- ts(reads.gr.pos)
            reads.gr.neg.ts <- ts(reads.gr.neg)
            ccf <- ccf(reads.gr.neg.ts, reads.gr.pos.ts, type = "correlation", plot=FALSE, lag.max=lag.max)
            keep <- ccf$lag>=0
            lag <- ccf$lag[keep]
            acf <- ccf$acf[keep]
            pos <- lag[which.max(acf)]
            if(plot){
                try({
                    plot(lag, acf, xlab="Lag", ylab="ACF", main=paste(f, "insertion size (not including reads length)"), type="l")
                    abline(v=pos, col="red")
                })
            }
            pos + mean(width(reads))
        }
    }, bamfiles, index, SIMPLIFY = TRUE)
    res
}