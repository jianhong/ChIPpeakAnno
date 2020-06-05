#' Output a summary of the occurrence of each pattern in the sequences.
#' 
#' Output a summary of the occurrence of each pattern in the sequences.
#' 
#' 
#' @param patternFilePath A character vector containing the path to the file to
#' read the patterns from.
#' @param format Either "fasta" (the default) or "fastq"
#' @param skip Single non-negative integer. The number of records of the
#' pattern file to skip before beginning to read in records.
#' @param BSgenomeName BSgenome object. Please refer to available.genomes in
#' BSgenome package for details
#' @param peaks \link[GenomicRanges:GRanges-class]{GRanges} containing the
#' peaks
#' @param outfile A character vector containing the path to the file to write
#' the summary output.
#' @param append TRUE or FALSE, default FALSE
#' @return A data frame with 3 columns as n.peaksWithPattern (number of peaks
#' with the pattern), n.totalPeaks (total number of peaks in the input) and
#' Pattern (the corresponding pattern). The summary will consider both strand
#' (including reverse complement).
#' @author Lihua Julie Zhu
#' @keywords misc
#' @export
#' @importFrom Biostrings readDNAStringSet reverseComplement
#' @importFrom utils write.table
#' @examples
#' 
#'     peaks = GRanges(seqnames=c("NC_008253", "NC_010468"), 
#'                     IRanges(start=c(100, 500), end=c(300, 600), 
#'                                names=c("peak1", "peak2")))
#'     filepath =system.file("extdata", "examplePattern.fa", 
#'                           package="ChIPpeakAnno")
#'     library(BSgenome.Ecoli.NCBI.20080805)
#'     summarizePatternInPeaks(patternFilePath=filepath, format="fasta", 
#'                             skip=0L, BSgenomeName=Ecoli, peaks=peaks)
#' 
summarizePatternInPeaks <-
  function (patternFilePath, format = "fasta", skip = 0L, BSgenomeName,
            peaks, outfile, append = FALSE)
  {
    if (missing(patternFilePath)) {
      stop("missing required parameter patternFilePath!")
    }
    if (!file.exists(patternFilePath)) {
      stop("patternFilePath specified as ", patternFilePath,
           " does not exsist!")
    }
    if (format != "fasta" && format != "fastq") {
      stop("format needs to be either fasta or fastq!")
    }
    if (missing(BSgenomeName) || !is(BSgenomeName, "BSgenome")) {
      stop("BSgenomeName is required as BSgenome object!")
    }
    if (missing(peaks) || (!is(peaks, "RangedData") && !is(peaks,
                                                           "GRanges"))) {
      stop("No valid peaks passed in. 
           It needs to be RangedData or GRanges object.")
    }
    if (is(peaks, "RangedData"))
      peaks <- toGRanges(peaks, format = "RangedData")
    if (!missing(outfile) && file.exists(outfile) && !append) {
      stop("outfile specified as ", outfile, 
           " already exists! Please rename the \n  
           outfile or set append = TRUE!")
    }
    seq = getAllPeakSequence(peaks, upstream = 0, downstream = 0,
                             genome = BSgenomeName)
    dict = readDNAStringSet(patternFilePath, format, use.names = TRUE)
    temp <- do.call(rbind, lapply(seq_along(dict), function(i) {
      getPosInSeqs(thispattern = dict[i], seq = seq, seqGRanges = peaks)
    }))
    if (!missing(outfile)) {
      write.table(temp, outfile, append = append, sep = "\t",
                  row.names = FALSE)
    }
    temp
  }


getPosInSeqs <-
  function(thispattern, seq, seqGRanges)
  {
    if (missing(thispattern) || !is(thispattern, "DNAStringSet")) {
      stop("thispattern is required as a DNAStringSet object!")
    }
    if (missing(seq)) {
      stop("No valid sequences passed in!")
    }
    patternName <- names(thispattern)
    revcomp.pattern = reverseComplement(thispattern)
    thispattern.t = as.character(thispattern)[[1]]
    thispattern.t = translatePattern(thispattern.t)
    revcomp.pattern.t = as.character(revcomp.pattern)[[1]]
    revcomp.pattern.t = translatePattern(revcomp.pattern.t)
    sequences <- seq$sequence
    seqInfo <- as.data.frame(seqGRanges)
    
    total <- do.call(rbind, lapply(seq_along(sequences), function(i) {
      pos.plus = gregexpr(thispattern.t, sequences[i], perl = TRUE)[[1]]
      pos.minus = gregexpr(revcomp.pattern.t, sequences[i], perl = TRUE)[[1]]
      if (pos.plus[1] > 0) {
        pattern.start <- as.numeric(pos.plus)
        pattern.end <- pattern.start + attr(pos.plus, "match.length") - 1
        patternSeq.found.plus <- 
          do.call(rbind, lapply(1:length(pattern.start), function(j){
          c(patternName, as.character(thispattern), 
            pattern.start[j], pattern.end[j], 
            substr(sequences[i], pattern.start[j], 
                   pattern.end[j]), "+", seqInfo[i,])
        }))
      }
      if (pos.minus[1] > 0) {
        pattern.start <- as.numeric(pos.minus)
        pattern.end <- pattern.start + attr(pos.minus, "match.length") - 1
        patternSeq.found.minus <- 
          do.call(rbind, lapply(1:length(pattern.start), function(j){
          c(patternName, as.character(thispattern), 
            pattern.start[j], pattern.end[j], 
            substr(sequences[i], pattern.start[j],
                   pattern.end[j]), "-",  seqInfo[i,])
        }))
        if (pos.plus[1] > 0)
          rbind(patternSeq.found.plus, patternSeq.found.minus)
        else
          patternSeq.found.minus
      }
      else if (pos.plus[1] > 0)
        patternSeq.found.plus
    }))
    if (length(total) == 0)
      cat( c(patternName, as.character(thispattern), 
             "not found in the input sequences!\n\n"))
    else
    {
      total <- as.data.frame(total)
      colnames(total)[1:6] <- c("motif name", "motif", "motifStartOffset",
                                "motifEndOffset", "motif found", 
                                "motifFoundStrand")
      chr <- as.character(total$seqnames)
      motifstart <- as.numeric(total$start) + as.numeric(total[,3]) -1
      motifend <- as.numeric(total$start) + as.numeric(total[,4]) -1
      total <- cbind(chr = chr, motifStart = motifstart, 
                     motifEnd = motifend, total)
      for (i in 1:dim(total)[2])
        total[,i] <- unlist(total[,i])
    }
    total
  }
