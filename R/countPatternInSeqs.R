countPatternInSeqs <- function(pattern, sequences){     
    if(missing(pattern) || class(pattern) != "DNAStringSet" )
    {
        stop("pattern is required as a DNAStringSet object!")
        
    }
    if (missing(sequences)) {
        stop("No valid sequences passed in!")
    }
    revcomp.pattern = reverseComplement(pattern)
    pattern = as.character(pattern)[[1]]
    pattern = translatePattern(pattern)
    revcomp.pattern = as.character(revcomp.pattern)[[1]]
    revcomp.pattern = translatePattern(revcomp.pattern)
    
    total =0 
    lapply(1:length(sequences), function(i){
        pos.plus = regexpr(pattern, sequences[i], perl=TRUE)[1]
        pos.minus = regexpr(revcomp.pattern, sequences[i], perl=TRUE)[1]
        if (pos.plus >0 || pos.minus >0)
        {
            total <<- total + 1;
        }
    })
    total
}