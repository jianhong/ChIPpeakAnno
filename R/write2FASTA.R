#' Write sequences to a file in fasta format
#' 
#' Write the sequences obtained from getAllPeakSequence to a file in fasta
#' format leveraging writeFASTA in Biostrings package. FASTA is a simple file
#' format for biological sequence data. A FASTA format file contains one or
#' more sequences and there is a header line which begins with a > proceeding
#' each sequence.
#' 
#' 
#' @param mySeq GRanges with varibles name and sequence ,e.g., results obtained
#' from getAllPeakSequence
#' @param file Either a character string naming a file or a connection open for
#' reading or writing. If "" (the default for write2FASTA), then the function
#' writes to the standard output connection (the console) unless redirected by
#' sink
#' @param width The maximum number of letters per line of sequence
#' @return Output as FASTA file format to the naming file or the console.
#' @author Lihua Julie Zhu
#' @keywords misc
#' @export
#' @importFrom Biostrings writeXStringSet
#' @examples
#' 
#' peaksWithSequences = GRanges(seqnames=c("1", "2"),
#' IRanges(start=c(1000, 2000), 
#' end=c(1010, 2010), 
#' names=c("id1", "id2")), 
#' sequence= c("CCCCCCCCGGGGG", "TTTTTTTAAAAAA"))
#' 
#' write2FASTA(peaksWithSequences, file="testseq.fasta", width=50)
#' 
write2FASTA <- function(mySeq, file="", width=80){
    if(!inherits(mySeq, c("GRanges"))){
        stop("mySeq must be an object of GRanges")
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
