#' Oligonucleotide enrichment analysis
#'
#' Test if the oligonucleotide in given region of chromosomes is enriched or not.
#'
#' @param filepath A character vector containing the path to the file to read the patterns from.
#'
#' @param format Either "fasta" (default) or "fastq".
#'
#' @param peaks GRanges containing the peaks.
#'
#' @param genomeName BSgenome object. Please refer to available.genomes in BSgenome package for details.
#'
#' @param background The method to get the background of compared oligonucleotide. "randomlyPick" (default) or "shuffle".
#'
#' @param chromosome Specify which chromosome will be selected to randomly pick back ground sequences. Default is the chromosome in peaks. Note that this parameter is valid for 'randomlyPick' method.
#'
#' @param k K-let size, dufault is 1, note that this parameter is valid for 'shuffle' method. See \link[universalmotif:shuffle_sequences]{k} for more information.
#'
#' @param rng.seed Set random number generator seed. See \link[universalmotif:shuffle_sequences]{rng.seed} for more information.
#'
#' @param method The method of shuffle sequence.'markov' (default),'euler','linear', only relevant if k> 1. See \link[universalmotif:shuffle_sequences]{method} for more information.
#'
#' @param nthreads nthreads = 0 uses all available threads. Note that no speed up will occur for jobs with only a single sequence, default is 1, note that this parameter is valid for 'shuffle' method. See \link[universalmotif:shuffle_sequences]{nthreads} for more information.
#'
#' @param window Logical value, Shuffle sequences iteratively over windows instead of all at once. Default is FALSE, note that this parameter is valid for 'shuffle' method. See \link[universalmotif:shuffle_sequences]{window} for more information.
#'
#' @param window.size Window size. Can be a fraction less than one, or an integer representing the actual window size. Default is 0.01, note that this parameter is valid for 'shuffle' method. See \link[universalmotif:shuffle_sequences]{window.size} for more information.
#'
#' @param window.overlap Overlap between windows. Can be a fraction less than one, or an integer representing the actual overlap size. Default is 0.01, note that this parameter is valid for 'shuffle' method. See \link[universalmotif:shuffle_sequences]{window.overlap} for more information.
#'
#' @param times The times of getting background sequence for the null distribution, default is 1000
#'
#' @param alpha The significant level, default is 0.05
#'
#' @return
#'
#' A data frame with 5 columns as x (number of match of the pattern), n (total number of oligonucleotide with the same length of pattern in the input) and prop.background (the proportions of pattern in background sequence), binom.pvalue (p value for the null that probabilities of success equal certain given values ) and threshold(p value threshold with given times of sample or shuffle).
#'
#' @details
#'
#' For the paramter k, method, nthreads, window, window.size and window.overlap, please see \link[universalmotif]{shuffle_sequences} for more information.
#'
#' @author Junhui Li
#'
#' @import universalmotif
#'
#' @examples
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' filepath =system.file("extdata", "examplePattern.fa", package="ChIPpeakAnno")
#' peaks = GRanges(seqnames=c("chr17","chr3","chr12","chr8"),
#'                  IRanges(start=c(41275784,10076141,4654135,31024288),
#'                          end=c(41276382,10076732,4654728,31024996),
#'                          names=paste0("peak",1:4)))
#'
#' result <- oligoNucleotideEnrichment(filepath=filepath,
#' peaks=peaks,
#' genomeName=Hsapiens,
#' background="randomPick")
#'
#' @export

oligoNucleotideEnrichment <- function(filepath,
                                      format="fasta",
                                      peaks,
                                      genomeName,
                                      background=c("randomPick","shuffle"),
                                      chromosome=NULL,
                                      k = 1,
                                      rng.seed = sample.int(10000, 1),
                                      method = c('markov','euler','linear'),
                                      nthreads = 1,
                                      window = FALSE,
                                      window.size = 0.1,
                                      window.overlap = 0.01,
                                      times=1000,
                                      alpha=0.05){
  background <- match.arg(background)
  method <- match.arg(method)
  actualSeqs <- getAllPeakSequence(peaks, upstream = 0, downstream = 0, genome = genomeName)
  actualP <- oligoNucleotideSummary(filepath,format=format,seqs=actualSeqs)
  nPeak <- length(peaks)
  allseqLen <- seqlengths(genomeName)
  seqWidth <- width(peaks)
  statDist <- matrix(NA,0,4)
  colnames(statDist) <- c("x","n","prop.background","binom.pvalue")
  if(background=="randomPick"){
    for(n in 1:times){
      if(is.null(chromosome)==TRUE){
        chrs <- as.character(seqnames(peaks))
        givenSeqLen <- allseqLen[as.character(seqnames(peaks))]
      }else{
        chrs <- sample(chromosome,nPeak)
        givenSeqLen <- allseqLen[chrs]
      }
      #i=1

      startPos <- NULL
      for(i in 1:nPeak){
        startPos <- append(startPos,sample(1:(givenSeqLen[i]-seqWidth[i]),1))
      }
      endPos <- startPos + seqWidth - 1
      backPeak <- GRanges(seqnames=chrs,
                          IRanges(start=startPos,
                                  end=endPos,
                                  names=paste0("peak",1:4)))
      backSeq <- getAllPeakSequence(backPeak, upstream = 0, downstream = 0, genome = genomeName)
      statTmp <- oligoNucleotideSummary(filepath,format=format,seqs=backSeq)
      rownames(statTmp) <- paste0(n,"_",rownames(statTmp))
      statDist <- rbind(statDist,statTmp)
    }
  }else if(background=="shuffle"){
    colnames(mcols(actualSeqs))[3] <- "orig.sequence"
    names(actualSeqs$orig.sequence) <- do.call(paste, c(as.data.frame(actualSeqs)[,1:3], sep="_"))
    inputSeq <- DNAStringSet(actualSeqs$orig.sequence)
    for(n in 1:times){
      sequences.shuffled <- universalmotif::shuffle_sequences(inputSeq, k = k,method = method,nthreads = nthreads,rng.seed = rng.seed,window = window, window.size = window.size,window.overlap = window.overlap)
      actualSeqs$sequence <- as.vector(sequences.shuffled)
      statTmp <- oligoNucleotideSummary(filepath,format=format,seqs=actualSeqs)
      rownames(statTmp) <- paste0(n,"_",rownames(statTmp))
      statDist <- rbind(statDist,statTmp)
    }
  }
  sortStatDist <- matrix(0,times,0)
  for(i in 1:(nrow(statDist)/times)){
    sortStatDist <- cbind(sortStatDist,data.frame(sort(as.numeric(statDist[seq(i,nrow(statDist),nrow(statDist)/times),4]))))
  }
  colnames(sortStatDist) <- 1:(nrow(statDist)/times)
  thresholdP <- sortStatDist[times*alpha,]
  actualP <- cbind(actualP,t(thresholdP))
  colnames(actualP)[5] <- "threshold"
  return(actualP)
}

oligoNucleotideSummary <- function(filepath,format="fasta",seqs){
  ## read fasta file and convert seq character and traslate pattern
  oligoNset = readDNAStringSet(filepath=filepath, format=format, use.names = TRUE)
  patternName <- names(oligoNset)
  revcomp.pattern = reverseComplement(oligoNset)
  thispattern = as.character(oligoNset)
  thispattern.t = translatePattern(thispattern)
  revcomp.pattern = as.character(revcomp.pattern)
  revcomp.pattern.t = translatePattern(revcomp.pattern)

  ## get frequency of single nuleotide in seqs
  singleNucCount <- colSums(oligonucleotideFrequency(DNAStringSet(seqs$sequence),width=1))
  singleNucFre <- singleNucCount/sum(singleNucCount)

  ## get expected frequency of fasta
  expFreList <- NULL
  for(i in 1:length(thispattern)){
    subPattern <- thispattern[i]
    forwardExpFre <- 1
    backwardExpFre <- 1
    #j=unlist(strsplit(subPattern,""))[1]
    for(j in unlist(strsplit(subPattern,""))){
      rcj <- reverseComplement(DNAString(j))
      if(j %in% c('A','T','C','G')){
        sbase <- singleNucFre[j]
        rbase <- singleNucFre[as.character(rcj)]
      }else if(j == 'N'){
        sbase <- 1
        rbase <- 1
      }else if(j == 'D'){
        sbase <- sum(singleNucFre[c("G","A","T")])
        rbase <- sum(singleNucFre[c("C","A","T")])
      }else if(j == 'V'){
        sbase <- sum(singleNucFre[c("G","A","C")])
        rbase <- sum(singleNucFre[c("G","C","T")])
      }else if(j == 'B'){
        sbase <- sum(singleNucFre[c("G","C","T")])
        rbase <- sum(singleNucFre[c("G","A","C")])
      }else if(j == 'H'){
        sbase <- sum(singleNucFre[c("C","A","T")])
        rbase <- sum(singleNucFre[c("G","A","T")])
      }else if(j == 'W'){
        sbase <- sum(singleNucFre[c("A","T")])
        rbase <- sum(singleNucFre[c("T","A")])
      }else if(j == 'S'){
        sbase <- sum(singleNucFre[c("G","C")])
        rbase <- sum(singleNucFre[c("G","C")])
      }else if(j == 'K'){
        sbase <- sum(singleNucFre[c("G","T")])
        rbase <- sum(singleNucFre[c("C","A")])
      }else if(j == 'M'){
        sbase <- sum(singleNucFre[c("C","A")])
        rbase <- sum(singleNucFre[c("G","T")])
      }else if(j == 'Y'){
        sbase <- sum(singleNucFre[c("C","T")])
        rbase <- sum(singleNucFre[c("G","A")])
      }else if(j == 'R'){
        sbase <- sum(singleNucFre[c("G","A")])
        rbase <- sum(singleNucFre[c("C","T")])
      }
      forwardExpFre <- forwardExpFre*sbase
      backwardExpFre <- backwardExpFre*rbase
    }
    names(forwardExpFre) <- names(backwardExpFre) <- names(thispattern[i])
    expFreList <- append(expFreList,forwardExpFre+backwardExpFre)
  }

  ## get real frequency of input seq in peak region
  XNmatrix <- matrix(0,length(patternName),3)
  rownames(XNmatrix) <- patternName
  colnames(XNmatrix) <- c("x","n","prop.background")
  uniqWidth <- unique(width(oligoNset))
  #i=6
  for(i in uniqWidth){
    oligoNVec <- colSums(oligonucleotideFrequency(DNAStringSet(seqs$sequence),width=i))
    for(j in which(width(oligoNset) %in% i)){
      XNmatrix[j,2] <- sum(oligoNVec)
      pos.plus = gregexpr(thispattern.t[j], names(oligoNVec), perl = TRUE)
      pos.minus = gregexpr(revcomp.pattern.t[j], names(oligoNVec), perl = TRUE)
      XNmatrix[j,1] <- sum(oligoNVec[which(sapply(pos.plus, "[[", 1) > 0)]) + sum(oligoNVec[which(sapply(pos.minus, "[[", 1) > 0)])
    }
  }

  XNmatrix[,3] <- expFreList
  binom.pvalue <- NULL
  for(i in 1:nrow(XNmatrix)){
    binom.result <- binom.test(XNmatrix[i,1],XNmatrix[i,2],XNmatrix[i,3],alternative="greater")
    binom.pvalue <- append(binom.pvalue,binom.result$p.value)
  }
  XNmatrix <- data.frame(XNmatrix,binom.pvalue)
  return(XNmatrix)
}
