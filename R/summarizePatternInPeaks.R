summarizePatternInPeaks <-
    function(patternFilePath, format="fasta", skip=0L, BSgenomeName, 
             peaks, outfile, append=FALSE)
    {
        if(missing(patternFilePath))
        {
            stop("missing required parameter patternFilePath!")
        }
        if (!file.exists(patternFilePath))
        {
            stop("patternFilePath specified as ", patternFilePath, 
                 " does not exsist!")
        }
        if (format != "fasta" && format != "fastq")
        {
            stop("format needs to be either fasta or fastq!")
        }
        if (missing(BSgenomeName) || class(BSgenomeName) != "BSgenome")
        {
            stop("BSgenomeName is required as BSgenome object!")
        }
        if (missing(peaks) || (class(peaks) != "RangedData" && 
                                   class(peaks) !="GRanges")) {
            stop("No valid peaks passed in. It needs to 
                 be RangedData or GRanges object.")
        }
        if(class(peaks)=="RangedData") 
            peaks <- toGRanges(peaks, format="RangedData")
        if (!missing(outfile) && file.exists(outfile) && !append)
        {
            stop("outfile specified as ", outfile, 
                 " already exists! Please rename the 
                 outfile or set append = TRUE!")
        }
        seq = getAllPeakSequence(peaks, upstream = 0, 
                                 downstream = 0, genome=BSgenomeName)
        n.peaks = length(peaks)
        
        dict = readDNAStringSet(patternFilePath, format, use.names=TRUE)
        
        temp = do.call(rbind, lapply(1:length(dict), function(i){
            total=countPatternInSeqs(pattern=dict[i], sequences=seq$sequence)
            c(total, n.peaks, as.character(unlist(dict[i])))
        }
        ))
        
        colnames(temp) = c("n.peaksWithPattern", "n.totalPeaks", "Pattern")
        if (!missing(outfile))
        {
            write.table(temp,outfile, append=append, sep="\t",row.names=FALSE)
        }
        temp
    }