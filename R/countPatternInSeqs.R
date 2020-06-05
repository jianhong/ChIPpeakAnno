#' Output total number of patterns found in the input sequences
#' 
#' Output total number of patterns found in the input sequences
#' 
#' 
#' @param pattern DNAstringSet object
#' @param sequences a vector of sequences
#' @return Total number of occurrence of the pattern in the sequences
#' @author Lihua Julie Zhu
#' @seealso summarizePatternInPeaks, translatePattern
#' @keywords misc
#' @export
#' @importFrom Biostrings reverseComplement
#' @examples
#'   library(Biostrings)
#'   filepath = 
#'     system.file("extdata", "examplePattern.fa", package="ChIPpeakAnno")
#'   dict = readDNAStringSet(filepath = filepath, format="fasta", 
#'                           use.names=TRUE)
#'   sequences = c("ACTGGGGGGGGCCTGGGCCCCCAAAT", 
#'                 "AAAAAACCCCTTTTGGCCATCCCGGGACGGGCCCAT", 
#'                 "ATCGAAAATTTCC")
#'   countPatternInSeqs(pattern=dict[1], sequences=sequences)
#'   countPatternInSeqs(pattern=dict[2], sequences=sequences)
#'   pattern = DNAStringSet("ATNGMAA")
#'   countPatternInSeqs(pattern=pattern, sequences=sequences)
#' 
countPatternInSeqs <- function(pattern, sequences){     
    if(missing(pattern) || !is(pattern, "DNAStringSet" ))
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
