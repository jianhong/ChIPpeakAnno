summarizePatternInPeaks <-
  
  function (patternFilePath, format = "fasta", skip = 0L, BSgenomeName, 
            peaks, outfile, append = FALSE) {
    
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
      
      stop("No valid peaks passed in. It needs to be RangedData or GRanges object.")
      
    }
    
    if (is(peaks, "RangedData")) 
      
      peaks <- toGRanges(peaks, format = "RangedData")
    
    if (!missing(outfile) && file.exists(outfile) && !append) {
      
      stop("outfile specified as ", outfile, 
           " already exists! Please rename the \n                 outfile or set append = TRUE!")
      
    }
    
    seq = getAllPeakSequence(peaks, upstream = 0, downstream = 0, 
                             
                             genome = BSgenomeName)
    
    n.peaks = length(peaks)
    
    name.peaks <- names(peaks)
    
    if (length(name.peaks) < n.peaks)
      
      name.peaks <- paste("peak", seq.int(n.peaks))
    
    dict = readDNAStringSet(patternFilePath, format, use.names = TRUE)
    
    temp = do.call(rbind, lapply(seq_along(dict), function(i) {
      
      total = getPosInSeqs(pattern = dict[i], sequences = seq$sequence)
      
      if (n.peaks == 1)
        
        c(as.character(total), name.peaks, as.character(unlist(dict[i])))
      
      else
        
        cbind(as.character(total), name.peaks, as.character(unlist(dict[i])))
      
    }))
    
    colnames(temp) = c("PatternPosition", "nameOfPeaks", 
                       
                       "Pattern")
    
    temp <- subset(temp, temp[,1] != "")
    
    if (!missing(outfile)) {
      
      write.table(temp, outfile, append = append, sep = "\t", 
                  
                  row.names = FALSE)
      
    }
    
    temp
    
  }



getPosInSeqs <-
  
  function (pattern, sequences){
    
    if (missing(pattern) || !is(pattern, "DNAStringSet")) {
      
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
    
    total = ""
    
    total <- lapply(seq_along(sequences), function(i) {
      
      pos.plus = gregexpr(pattern, sequences[i], perl = TRUE)[[1]]
      
      pos.minus = gregexpr(revcomp.pattern, sequences[i], perl = TRUE)[[1]]
      
      if (pos.plus[1] > 0) {
        total = paste(total,  paste(as.character(pos.plus), ","), sep="plus:")
      }
      
      
      
      if ( pos.minus[1] > 0) {
        total = paste(total,  paste(as.character(pos.minus), ","), sep="minus:")
      }
      
      total
      
    })
    
    total
    
  }
