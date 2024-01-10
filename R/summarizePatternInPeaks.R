#' Output a summary of the occurrence and enrichment of each pattern in the 
#' sequences.
#' 
#' Output a summary of the occurrence and enrichment of each pattern in the 
#' sequences.
#'
#' @param patternFilePath Character value. The path to the file that contains 
#' the pattern. 
#' 
#' @param format Character value. The format of file containing the 
#' oligonucleotide pattern, either "fasta" (default) or "fastq".
#' 
#' @param BSgenomeName Character value. BSgenome object. Please refer to 
#' available.genomes in BSgenome package for details.
#' 
#' @param peaks Character value. \link[GenomicRanges:GRanges-class]{GRanges} 
#' containing the peaks.
#' 
#' @param revcomp Boolean value, if TURE, also search the reverse compliment of 
#' pattern. Default is TRUE.
#' 
#' @param method Character value. Method for pattern enrichment test, 
#' 'binom.test' (default) or 'permutation.test'.
#' 
#' @param expectFrequencyMethod Character value. Method for calculating the 
#' expected probability of pattern occurrence, 'Markov' (default) or 'Naive'.
#' 
#' @param MarkovOrder Integer value. The order of Markov chain. Default is 3.
#' 
#' @param bgdForPerm Character value. The method for obtaining the background 
#' sequence. 'chromosome' (default) selects background chromosome 
#' from chromosomes, refer to 'chromosome' parameter; 'shuffle' will obtain the
#' backgroud sequence by shufflubg any k-mers in peak sequences, refer to '...'.
#'
#' @param chromosome Character value. Relevant if "bgdForPerm='chromosome'".
#' 'asPeak' means to use the same chromosomes in peaks; 'random' means to use 
#' all chromosomes randomly. Default is 'asPeak'. 
#' 
#' @param nperm Integer value. The number of permutation test, default is 1000.
#'
#' @param alpha Numeric value. The significant level for permutation test, 
#' default is 0.05.
#'
#' @param ... Aditional parameter passed to function \\link[universalmotif]
#' {shuffle_sequences}
#'
#' @return A list including two data frames named 'motif_enrichment' and 
#' 'motif_occurrence'. The 'motif_enrichment' has four columns: 
#' \itemize{
#'  \item{"patternNum"}: {number of matched pattern}
#'  \item{"totalNumPatternWithSameLen"}:  {total number of pattern with the 
#'  same length}
#'  \item{"expectedRate"}:  {expected rate of pattern for 'binom.test' method}
#'  \item{"patternRate"}:  {real rate of pattern for 'permutation.test' method}
#'  \item{"pValueBinomTest"}:  {p value of bimom test for 'binom.test' method}
#'  \item{"cutOffPermutationTest"}:  {cut off of permutation test for 
#'  'permutation.test' method}
#' }
#' 
#'The 'motif_occurrence' has 14 columns:
#'
#' \itemize{
#'  \item{"motifChr"}: {Chromosome of motif}
#'  \item{"motifStartInChr"}:  {motif start position in chromosome}
#'  \item{"motifEndInChr"}:  {motif end position in chromosome}
#'  \item{"motifName"}:  {motif name}
#'  \item{"motifPattern"}:  {motif pattern}
#'  \item{"motifStartInPeak"}:  {motif start position in peak}
#'  \item{"motifEndInPeak"}: {motif end position in peak}
#'  \item{"motifFound"}:  {specific motif Found in peak}
#'  \item{"motifFoundStrand"}:  {strand of specific motif Found in peak, "-" 
#'  means reverse complement of motif found in peaks}
#'  \item{"peakChr"}:  {Chromosome of peak}
#'  \item{"peakStart"}:  {peak start position}
#'  \item{"peakEnd"}:  {peak end position}
#'  \item{"peakWidth"}:  {peak width}
#'  \item{"peakStrand"}:  {peak strand}
#' }
#' 
#' @details Please see \link[universalmotif]{shuffle_sequences} for the more 
#' information bout 'shuffle' method.
#' 
#' @author Lihua Julie Zhu, Junhui Li, Kai Hu
#' 
#' @export
#' 
#' @importFrom Biostrings readDNAStringSet reverseComplement
#' 
#' @importFrom Biostrings oligonucleotideFrequency
#' 
#' @importFrom stats binom.test
#' 
#' @importFrom universalmotif shuffle_sequences
#' 
#' @importFrom data.table data.table
#' 
#' @importFrom stringr str_split
#' 
#' @importFrom dplyr %>% pull left_join
#' 
#' @importFrom tidyr separate_rows
#' 
#' @importFrom tibble tibble
#' 
#' @keywords misc
#' 
#' @examples
#'                             
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' filepath <- system.file("extdata", "examplePattern.fa", 
#'                         package = "ChIPpeakAnno")
#' peaks <- GRanges(seqnames = c("chr17", "chr3", "chr12", "chr8"),
#'                  IRanges(start = c(41275784, 10076141, 4654135, 31024288),
#'                          end = c(41276382, 10076732, 4654728, 31024996),
#'                          names = paste0("peak", 1:4)))
#' result <- summarizePatternInPeaks(patternFilePath = filepath, peaks = peaks,
#'                                   BSgenomeName = Hsapiens)
#' 
summarizePatternInPeaks <- function(patternFilePath, 
                                    format = "fasta", 
                                    BSgenomeName,
                                    peaks,
                                    revcomp = TRUE,
                                    method = c("binom.test","permutation.test"),
                                    expectFrequencyMethod = c("Markov","Naive"),
                                    MarkovOrder = 3L,
                                    bgdForPerm = c("shuffle","chromosome"),
                                    chromosome = c("asPeak","random"),
                                    nperm = 1000,
                                    alpha = 0.05,
                                    ...) {
  method <- match.arg(method)
  expectFrequencyMethod <- match.arg(expectFrequencyMethod)
  bgdForPerm <- match.arg(bgdForPerm)
  chromosome <- match.arg(chromosome)
  
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

  seqs <- getAllPeakSequence(peaks, upstream = 0, downstream = 0, 
                             genome = BSgenomeName)
  
  # seqs <- getAllPeakSequence(peaks, BSgenome.Hsapiens.UCSC.hg19, upstream = 0, downstream = 0)
  # patternVec <- DNAStringSet("AACCCA")
  # names(patternVec) <- "testRun"
  # revcomp <- FALSE
  
  patternVec <- readDNAStringSet(patternFilePath, format, use.names = TRUE)
  patternName <- names(patternVec)
  temp <- do.call(rbind, lapply(seq_along(patternVec), function(i) {
    getPosInSeqs(thispattern = patternVec[i], seq = seqs, revcomp = revcomp)
  }))
  temp2 <- patternEnrichmentInPeak(patternVec,
                                   revcomp = revcomp,
                                   BSgenomeName = BSgenomeName,
                                   seqs = seqs, 
                                   method = method, 
                                   expectFrequencyMethod = expectFrequencyMethod,
                                   bgdForPerm = bgdForPerm,
                                   chromosome = chromosome,
                                   MarkovOrder = MarkovOrder,
                                   nperm = nperm,
                                   alpha = alpha,
                                   ...)
  
  output <- list()
  output["motif_enrichment"] <- list(temp2)
  output["motif_occurrence"] <- list(temp)
  return(output)
}

patternEnrichmentInPeak <- function(patternVec,
                                    revcomp = TRUE,
                                    seqs,
                                    method = c("binom.test", "permutation.test"),
                                    expectFrequencyMethod = c("Markov", "Naive"),
                                    bgdForPerm = c("chromosome", "shuffle"),
                                    chromosome = c("asPeak","random"),
                                    BSgenomeName,
                                    MarkovOrder = 3L,
                                    nperm = 1000,
                                    alpha = 0.05,
                                    ...) {
  method <- match.arg(method)
  expectFrequencyMethod <- match.arg(expectFrequencyMethod)
  bgdForPerm <- match.arg(bgdForPerm)
  stopifnot("The 'nperm' or 'alpha' parameter should be increased to a 
            sufficient extent."=round(nperm*alpha,0) != 0)
  
  pattern <- as.character(patternVec)
  if(method == "binom.test"){
    peakPvalue <- binomEnrichment(pattern, 
                                  expectFrequencyMethod = expectFrequencyMethod,
                                  revcomp = revcomp, 
                                  seqs = seqs, 
                                  MarkovOrder = MarkovOrder)
  }else if(method == "permutation.test"){
    peakPvalue <- permutationEnrichment(pattern, 
                                        revcomp = revcomp, 
                                        seqs = seqs, 
                                        BSgenomeName = BSgenomeName, 
                                        bgdForPerm = bgdForPerm,
                                        chromosome = chromosome,
                                        nperm = nperm, 
                                        alpha = alpha,
                                        ...)
  }
  return(peakPvalue)
}

binomEnrichment <- function(pattern, 
                            expectFrequencyMethod = c("Naive","Markov"), 
                            revcomp = TRUE, 
                            seqs, 
                            MarkovOrder = 3L) {
  ## get frequency of single nucleotide in seqs
  expFreqSet <- expectMotifFrequency(pattern, 
                                     revcomp = revcomp, 
                                     seqs = seqs, 
                                     expectFrequencyMethod = expectFrequencyMethod,
                                     MarkovOrder = MarkovOrder)
  ## get real frequency of input seq in peak region
  peakFreq <- matrix(0,length(pattern),3)
  rownames(peakFreq) <- pattern
  colnames(peakFreq) <- c("patternNum","totalNumPatternWithSameLen",
                          "expectedRate")
  peakFreq[,3] <- expFreqSet[,1]
  
  ## get expected frequency of fasta
  pattern.t <- translatePattern(pattern)
  backwardpattern.t <- NULL
  backwardpattern <- NULL
  if(revcomp == TRUE){
    backwardpattern <- reverseComplement(DNAStringSet(pattern))
    backwardpattern <- as.character(backwardpattern)
    backwardpattern.t <- translatePattern(backwardpattern)
  }
  peakFreq[,1:2] <- t(sapply(seq.int(pattern),function(i){
    oligoNVec <- colSums(oligonucleotideFrequency(DNAStringSet(seqs$sequence),
                                                  width=nchar(pattern[i])))
    posPlus = gregexpr(pattern.t[i], names(oligoNVec), perl = TRUE)
    posPlusFreq <- sum(oligoNVec[which(sapply(posPlus, "[[", 1) > 0)])
    posMinusFreq <- 0
    if(revcomp == TRUE){
      posMinus = gregexpr(backwardpattern.t[i], names(oligoNVec), perl = TRUE)
      posMinusFreq <- sum(oligoNVec[which(sapply(posMinus, "[[", 1) > 0)])
    }
    c(posPlusFreq + posMinusFreq, sum(oligoNVec))
  }))
  
  pValueBinomTest <- sapply(seq.int(nrow(peakFreq)),function(i){
    stats::binom.test(peakFreq[i,1],peakFreq[i,2],peakFreq[i,3],
                      alternative="greater")$p.value
  })
  peakFreq <- data.frame(peakFreq,pValueBinomTest)
  return(peakFreq)
}

expectMotifFrequency <- function(pattern, 
                                 revcomp = TRUE, 
                                 seqs,
                                 expectFrequencyMethod = c("Naive","Markov"),
                                 MarkovOrder = 3L) {
  
  if(expectFrequencyMethod == "Naive"){
    ACGTcount <- colSums(oligonucleotideFrequency(DNAStringSet(seqs$sequence),
                                                  width=1))
    ACGTfreq <- ACGTcount/sum(ACGTcount)
    iuapc <- data.table(code = c("A", "C", "G", "T", "R", "Y", "S", "W", "K", 
                                 "M", "B", "D", "H", "V", "N"),
                        base = c("A", "C", "G", "T", "AG", "CT", "GC", "AT", 
                                "GT", "AC", "CGT", "AGT", "ACT", "ACG", "ACGT"))
    allBaseFreq <- sapply(1:nrow(iuapc),function(i){
      sum(ACGTfreq[seqinr::s2c(as.character(iuapc[i,2]))])
    })
    names(allBaseFreq) <- as.matrix(iuapc)[,1]
    
    ## get expected frequency of forward and backward fasta 
    pattern.t <- translatePattern(pattern)
    backwardpattern.t <- NULL
    backwardpattern <- NULL
    if(revcomp == TRUE){
      backwardpattern <- reverseComplement(DNAStringSet(pattern))
      backwardpattern <- as.character(backwardpattern)
      backwardpattern.t <- translatePattern(backwardpattern)
    }
    expFreqSet <- sapply(c(pattern,backwardpattern),function(i){
      patternTab <- table(strsplit(i,"")[[1]])
      expFre <- prod(allBaseFreq[names(patternTab)]^patternTab)
      expFre
    })
    expFreqSet <- colSums(matrix(expFreqSet,ncol=length(pattern),byrow = TRUE))
    names(expFreqSet) <- pattern
    expectedFreq <- data.frame(expFreqSet)
    colnames(expectedFreq) <- "expected_frequency"
  }else if(expectFrequencyMethod == "Markov"){
    
    oligoLength <- nchar(pattern)
    stopifnot("motif length should > 3 and < 13 with Markov method"=
                all(oligoLength>3&oligoLength<13))
    MarkovOrder <- as.integer(MarkovOrder)
    stopifnot("MarkovOrder should > 0 and < 6" = MarkovOrder > 0 & MarkovOrder < 6)
    stopifnot("MarkovOrder should be less than motif length - 2"=
                all(MarkovOrder < oligoLength - 2))
    # stopifnot("The seqinr package is required." = 
    #             requireNamespace("seqinr", quietly = TRUE))
    if(is(seqs, "GRanges")){
      seqs_content <- seqs$sequence
    }else{
      if(is(seqs, "DNAStringSet")){
        seqs_content <- as.character(seqs)
      }else{
        if(!is.character(seqs)){
          stop("seqs must be an object of DNAStringSet or DNAString",
               "or output of getAllPeakSequence")
        }
      }
    }
    sequence <- tolower(seqs_content)
    len <- length(sequence)
    #i=6
    expectedFreq <- sapply(unique(oligoLength),function(i){
      oligoWords <- seqinr::words(i)
      dict <- PDict(DNAStringSet(oligoWords))
      sequence.tbl <- 1
      maxsize <- 100000
      if(len > maxsize){
        sequence.tbl <- table(sequence)
        if(length(sequence.tbl)>maxsize){
          stop("Can not handle such huge dataset.")
        }
        cnt <- vcountPDict(dict, subject=DNAStringSet(names(sequence.tbl)), 
                           max.mismatch=0, min.mismatch=0, 
                           with.indels=FALSE, fixed=TRUE)
      }else{
        cnt <- vcountPDict(dict, subject=DNAStringSet(sequence), 
                           max.mismatch=0, min.mismatch=0, 
                           with.indels=FALSE, fixed=TRUE)
      }
      
      rownames(cnt) <- oligoWords
      cnt <- t(cnt)
      mergeRevcomp <- function(mat){
        coln <- colnames(mat)
        index <- nrow(mat)
        map <- c(a="t", c="g", g="c", t="a")
        revComp <- function(.ele){
          paste(map[rev(seqinr::s2c(.ele))], collapse="")
        }
        coln.rev <- sapply(coln, revComp)
        coln.id <- 1:length(coln)
        coln.rev.id <- match(coln.rev, coln)
        coln.ids <- apply(cbind(coln.id, coln.rev.id), 1, sort)
        coln.ids <- unique(t(coln.ids))
        mat1 <- mat[, coln.ids[, 1]] > 0
        mat2 <- mat[, coln.ids[, 2]] > 0
        if (index == 1) {
          mat1 <- t(as.data.frame(mat1))
          mat2 <- t(as.data.frame(mat2))
        }
        mat1 <- colSums(mat1)
        mat2 <- colSums(mat2)
        coln <- ifelse(mat1>mat2, names(mat1), names(mat2))
        mat <- mat[, coln.ids[, 1]] + mat[, coln.ids[, 2]]
        if (index == 1) {
          mat <- t(as.data.frame(mat))
        }
        colnames(mat) <- coln
        mat
      }
      
      if(revcomp) cnt <- mergeRevcomp(cnt)
      
      #cnt <- cnt>0
      mode(cnt) <- "logical"
      mode(cnt) <- "numeric"
      if(len > maxsize) cnt <- cnt * as.numeric(sequence.tbl)
      cntSum <- colSums(cnt)
      freqs <- oligoFrequency(sequence, MarkovOrder=MarkovOrder)
      namesFreqs <- unique(c(seqinr::words(1), 
                             seqinr::words(MarkovOrder),
                             seqinr::words(MarkovOrder+1)))
      
      names(freqs) <- tolower(names(freqs))
      stopifnot(all(namesFreqs %in% names(freqs)))
      f <- sapply(names(cntSum), function(.ele) {
        m1 <- substring(.ele, 1:(i-MarkovOrder),
                        (MarkovOrder+1):i)
        m0 <- substring(.ele, 2:(i-MarkovOrder),
                        (MarkovOrder+1):(i-1))
        prod(freqs[m1])/prod(freqs[m0])
      })
      #i=6
      allPatternList <- expandPattern(pattern[oligoLength %in% i])
      
      expFreqSubset <- sapply(allPatternList,function(seq){
        sum(f[names(f) %in% tolower(seq)])
      })
      expFreqSubset
    })
    expectedFreq <- as.data.frame(expectedFreq)
    colnames(expectedFreq) <- "expected_frequency"
  }
  return(expectedFreq)
}


expandPattern <- function(pattern){
  pattern = toupper(pattern)
  iuapc <- data.table(code = c("A", "C", "G", "T", "R", "Y", "S", "W", "K", "M",
                               "B", "D", "H", "V", "N"),
                      base = c("A", "C", "G", "T", "AG", "CT", "GC", "AT", "GT",
                               "AC", "CGT", "AGT", "ACT", "ACG", "ACGT"))
  
  allExpPattern <- lapply(pattern,function(seq){
    allSeqMat <- tibble(seq) %>%
      separate_rows(seq, sep = '(?<=.)(?=.)') %>%
      left_join(iuapc, by = c("seq" = "code")) %>%
      pull(base) %>%
      str_split("") %>%
      expand.grid(stringsAsFactors = FALSE)
    apply(allSeqMat,1,seqinr::c2s)
  })
  names(allExpPattern) <- pattern
  allExpPattern
}

permutationEnrichment <- function(pattern, 
                                  revcomp = TRUE, 
                                  seqs, 
                                  BSgenomeName, 
                                  bgdForPerm = c("chromosome","shuffle"), 
                                  chromosome = c("asPeak","random"),
                                  nperm = 1000, 
                                  alpha = 0.05,
                                  ...) {
  pattern.t <- translatePattern(pattern)
  backwardpattern.t <- NULL
  backwardpattern <- NULL
  if (revcomp == TRUE) {
    backwardpattern <- reverseComplement(DNAStringSet(pattern))
    backwardpattern <- as.character(backwardpattern)
    backwardpattern.t <- translatePattern(backwardpattern)
  }
  allseqLen <- seqlengths(BSgenomeName)
  seqWidth <- width(seqs)
  nPeak <- length(seqs)
  colnames(mcols(seqs))[4] <- "orig.sequence"
  names(seqs$orig.sequence) <- do.call(paste, c(as.data.frame(seqs)[, 1:3], 
                                                sep = "_"))
  inputSeq <- DNAStringSet(seqs$orig.sequence)
  candiChrom <- names(allseqLen >= max(seqWidth))
  stopifnot("The length of peak sequence should be less than the length of 
  chromosomes"=all(seqWidth < allseqLen[as.character(seqnames(seqs)@values)]))
  
  peakFreq <- matrix(0,length(pattern),4)
  rownames(peakFreq) <- pattern
  colnames(peakFreq) <- c("patternNum","totalNumPatternWithSameLen",
                          "patternRate","cutOffPermutationTest")
  peakFreq[,1:2] <- t(sapply(seq.int(pattern),function(i){
    oligoNVec <- colSums(oligonucleotideFrequency(
      DNAStringSet(seqs$orig.sequence),width=nchar(pattern[i])))
    posPlus = gregexpr(pattern.t[i], names(oligoNVec), perl = TRUE)
    posPlusFreq <- sum(oligoNVec[which(sapply(posPlus, "[[", 1) > 0)])
    posMinusFreq <- 0
    if(revcomp == TRUE){
      posMinus = gregexpr(backwardpattern.t[i], names(oligoNVec), perl = TRUE)
      posMinusFreq <- sum(oligoNVec[which(sapply(posMinus, "[[", 1) > 0)])
    }
    c(posPlusFreq + posMinusFreq, sum(oligoNVec))
  }))
  peakFreq[,3] <- peakFreq[,1]/peakFreq[,2]
  
  ep <- list(...)
  permutationFreq <- do.call(rbind, lapply(seq.int(nperm), 
                                           function(n, ep) {
                                             if (bgdForPerm == "chromosome") {
      if (chromosome=="random"){
        chrs <- sample(candiChrom, nPeak, replace = TRUE)
        givenSeqLen <- allseqLen[chrs]
      } else if(chromosome=="asPeak"){
        chrs <- as.character(seqnames(seqs))
        givenSeqLen <- allseqLen[chrs]
      }
      startPos <- vapply(givenSeqLen-seqWidth,function(x){sample(seq.int(x),1)},
                         numeric(1))
      endPos <- startPos + seqWidth - 1
      backgroudPeak <- GRanges(seqnames = chrs,
                               IRanges(start = startPos,
                                       end = endPos,
                                       names = paste0("peak",seq.int(nPeak))))
      backgroudPeakseq <- getAllPeakSequence(backgroudPeak, upstream = 0,
                                             downstream = 0, 
                                             genome = BSgenomeName)
    
      } 
                                             else if (bgdForPerm == "shuffle") {
      sequences.shuffled <- do.call(universalmotif::shuffle_sequences, c(list(inputSeq), ep))
      seqs$sequence <- as.vector(sequences.shuffled)
      backgroudPeakseq <- seqs
    }
  ## get frequency from permutation
  sapply(seq.int(pattern),function(i){
  oligoNVec <- colSums(oligonucleotideFrequency(DNAStringSet(
    backgroudPeakseq$sequence),width=nchar(pattern[i])))
  posPlus <- gregexpr(pattern.t[i], names(oligoNVec), perl = TRUE)
  expFre <- sum(oligoNVec[which(sapply(posPlus, "[[", 1) > 0)])
  if(revcomp == TRUE){
    posMinus <- gregexpr(backwardpattern.t[i], names(oligoNVec), perl = TRUE)
    expFreMinus <- sum(oligoNVec[which(sapply(posMinus, "[[", 1) > 0)])
    expFre <- expFre + expFreMinus
  }
  c(expFre/sum(oligoNVec))
})
                                             }, ep = ep))
  colnames(permutationFreq) <- names(pattern)
  permutationFreqSorted <- apply(permutationFreq, 2, sort, decreasing=TRUE)
  peakFreq[,4] <- as.vector(permutationFreqSorted[round(nperm*alpha,0),])
  return(peakFreq)
}

getPosInSeqs <- function(thispattern, seq, revcomp = TRUE){
  if (missing(thispattern) || !is(thispattern, "DNAStringSet")) {
    stop("thispattern is required as a DNAStringSet object!")
  }
  if (missing(seq)) {
    stop("No valid sequences passed in!")
  }
  patternName <- names(thispattern)
  thispattern.t = as.character(thispattern)[[1]]
  thispattern.t = translatePattern(thispattern.t)
  if(revcomp == TRUE){
    revcomp.pattern = reverseComplement(thispattern)
    revcomp.pattern.t = as.character(revcomp.pattern)[[1]]
    revcomp.pattern.t = translatePattern(revcomp.pattern.t)
  }
  sequences <- seq$sequence
  #seqInfo <- as.data.frame(seqGRanges)
  seqInfo <- as.data.frame(seq)[,1:5]
  
  total <- do.call(rbind, lapply(seq_along(sequences), function(i) {
    pos.plus = gregexpr(thispattern.t, sequences[i], perl = TRUE)[[1]]
    if (pos.plus[1] > 0) {
      pattern.start <- as.numeric(pos.plus)
      pattern.end <- pattern.start + attr(pos.plus, "match.length") - 1
      patternSeq.found.plus <- 
        do.call(rbind, lapply(1:length(pattern.start), function(j){
          c(patternName, as.character(thispattern), 
            pattern.start[j], pattern.end[j], 
            substr(sequences[i], pattern.start[j], 
                   pattern.end[j]),"+", seqInfo[i,])
        }))
    }
    if(revcomp == TRUE){
      pos.minus = gregexpr(revcomp.pattern.t, sequences[i], perl = TRUE)[[1]]
      if(pos.minus[1] > 0) {
        pattern.start <- as.numeric(pos.minus)
        pattern.end <- pattern.start + attr(pos.minus, "match.length") - 1
        patternSeq.found.minus <- 
          do.call(rbind, lapply(1:length(pattern.start), function(j){
            c(patternName, as.character(thispattern), 
              pattern.start[j], pattern.end[j], 
              substr(sequences[i], pattern.start[j],
                     pattern.end[j]), "-",  seqInfo[i,])
          }))
        if(pos.plus[1] > 0){
          rbind(patternSeq.found.plus, patternSeq.found.minus)
        }else{
          patternSeq.found.minus
        }
      }else{
        if(pos.plus[1] > 0){
          patternSeq.found.plus
        }
      }
    }else{
      if(pos.plus[1] > 0){
        patternSeq.found.plus
      }
    }
  }))
  if (length(total) == 0){
    cat( c(patternName, as.character(thispattern), 
           "not found in the input sequences!\n\n"))
  }else{
    total <- as.data.frame(total,stringsAsFactors = FALSE)
    colnames(total)[1:6] <- c("motifName", "motifPattern", "motifStartInPeak",
                              "motifEndInPeak", "motifFound", 
                              "motifFoundStrand")
    motifstartInChromosome <- as.numeric(total$start) + as.numeric(total[,3]) -1
    motifendInChromosome <- as.numeric(total$start) + as.numeric(total[,4]) -1
    total <- cbind(motifChr = motifstartInChromosome,
                   motifStartInChr = motifstartInChromosome, 
                   motifEndInChr = motifendInChromosome, 
                   total)
    colnames(total)[(ncol(total)-4):ncol(total)] <- c("peakChr",
                                                      "peakStart",
                                                      "peakEnd",
                                                      "peakWidth",
                                                      "peakStrand")
    for (i in 1:dim(total)[2]) total[,i] <- unlist(total[,i])
    total$motifChr <- total$peakChr
  }
  total
}


