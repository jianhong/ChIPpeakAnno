write2FASTA <- function(mySeq, file="", width=80){
    if(!inherits(mySeq, c("RangedData", "GRanges"))){
        stop("mySeq must be an object of RangedData or GRanges")
    }
    if(is(mySeq, "RangedData")){
        mySeq <- toGRanges(mySeq)
    }
    if(is.null(mySeq$sequence)){
        stop("metadata must contain sequence")
    }
    if(!is.null(names(mySeq))){
        descriptions <- names(mySeq)
    }else{
        descriptions <- 
            paste("X", 
                  formatC(1:length(mySeq), 
                          width=nchar(as.character(length(mySeq))),
                          flag="0"), 
                  as.character(seqnames(mySeq)), ":", 
                  start(mySeq), "-",
                  end(mySeq), ":",
                  as.character(strand(mySeq)))
    }
    
    sequences <- mySeq$sequence
    names(sequences) <- descriptions
    writeXStringSet(as(sequences,"XStringSet"), file, width=width)
}
