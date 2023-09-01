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
#' @param genome BSgenome object or mart object. Please refer to available.genomes in BSgenome package and useMart in bioMaRt package for details
#'
#' @param background The method to get the background of compared oligonucleotide. "randomlyPick" (default) or "shuffle".
#'
#' @param chromosome Specify which chromosome will be selected to randomly pick back ground sequences. Default is the chromosome in peaks. Note that this parameter is valid for 'randomlyPick' method.
#'
#' @param ... could be parameters of function \code{\link[universalmotif]{shuffle_sequences}}
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
#' @import  Biostrings
#' 
#' @importFrom stats binom.test
#' 
#' @importFrom universalmotif shuffle_sequences
#'
#' @examples
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' filepath =system.file("extdata", "examplePattern.fa", package="ChIPpeakAnno")
#' peaks = GRanges(seqnames=c("chr17","chr3","chr12","chr8"),
#'                  IRanges(start=c(41275784,10076141,4654135,31024288),
#'                          end=c(41276382,10076732,4654728,31024996),
#'                          names=paste0("peak",1:4)))
#' \dontrun{
#' result <- oligoNucleotideEnrichment(filepath=filepath,
#' peaks=peaks,
#' genome=Hsapiens,
#' background="randomPick")
#' }
#' @export

oligoNucleotideEnrichment <- function(filepath,
                                      format="fasta",
                                      peaks,
                                      genome,
                                      background=c("randomPick","shuffle"),
                                      chromosome=NULL,
                                      ...,
                                      times=1000,
                                      alpha=0.05){
  if(!file.exists(filepath)){
    stop("file doesn't exist!")
  }
  background <- match.arg(background)
  method <- match.arg(method)
  if(!is(genome, "BSgenome") & !is(genome, "Mart") | missing(genome)){
    stop("genome is required parameter, 
             please pass in either a BSgenome object or a Mart object!")
  }
  if(round(times*alpha,0)==0){
    stop("The paramter 'times' should be increased to large enough!")
  }
  nPeak <- length(peaks)
  allseqLen <- seqlengths(genome)
  seqWidth <- width(peaks)
  candiChrom <- names(allseqLen >= max(seqWidth))
  if(any(seqWidth > allseqLen[as.character(seqnames(peaks)@values)])){
    stop("The length of sequence should be less than the length of chromosomes")
  }
  actualSeqs <- getAllPeakSequence(peaks, upstream = 0, downstream = 0, genome = genome)
  actualP <- oligoNucleotideSummary(filepath,format=format,seqs=actualSeqs)
  
  statDist <- matrix(NA,0,4)
  colnames(statDist) <- c("x","n","prop.background","binom.pvalue")
  if(background=="randomPick"){
    statDist <- do.call(rbind,lapply(seq.int(times),function(n) {
      if(is.null(chromosome)==TRUE){
        chrs <- as.character(seqnames(peaks))
        givenSeqLen <- allseqLen[as.character(seqnames(peaks))]
      }else{
        chrs <- sample(candiChrom,nPeak)
        givenSeqLen <- allseqLen[chrs]
      }
      startPos <- sapply(givenSeqLen-seqWidth,function(x){sample(seq.int(x),1)})
      endPos <- startPos + seqWidth - 1
      backPeak <- GRanges(seqnames=chrs,
                          IRanges(start=startPos,
                                  end=endPos,
                                  names=paste0("peak",seq.int(nPeak))))
      backSeq <- getAllPeakSequence(backPeak, upstream = 0, downstream = 0, genome = genome)
      statTmp <- oligoNucleotideSummary(filepath,format=format,seqs=backSeq)
      rownames(statTmp) <- paste0(n,"_",rownames(statTmp))
      statTmp
    }))
  }else if(background=="shuffle"){
    colnames(mcols(actualSeqs))[3] <- "orig.sequence"
    names(actualSeqs$orig.sequence) <- do.call(paste, c(as.data.frame(actualSeqs)[,1:3], sep="_"))
    inputSeq <- DNAStringSet(actualSeqs$orig.sequence)
    statDist <- do.call(rbind,lapply(seq.int(times),function(n){
      sequences.shuffled <- universalmotif::shuffle_sequences(inputSeq, ...)
      actualSeqs$sequence <- as.vector(sequences.shuffled)
      statTmp <- oligoNucleotideSummary(filepath,format=format,seqs=actualSeqs)
      rownames(statTmp) <- paste0(n,"_",rownames(statTmp))
      statTmp
    }))
  }
  sortStatDist <- do.call(cbind,lapply(seq.int(nrow(actualP)),function(i){
    sort(as.numeric(statDist[seq(i,nrow(statDist),nrow(actualP)),4]))
  }))
  colnames(sortStatDist) <- seq.int(nrow(actualP))
  thresholdP <- sortStatDist[round(times*alpha,0),]
  actualP <- cbind(actualP,thresholdP)
  colnames(actualP)[5] <- "threshold"
  return(actualP)
}

oligoNucleotideSummary <- function(filepath,format="fasta",seqs){
  if(!file.exists(filepath)){
    stop("file doesn't exist!")
  }
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
  D <- sum(singleNucFre[c("G","A","T")])
  H <- sum(singleNucFre[c("C","A","T")])
  B <- sum(singleNucFre[c("G","C","T")])
  V <- sum(singleNucFre[c("G","A","C")])
  W <- sum(singleNucFre[c("A","T")])
  S <- sum(singleNucFre[c("G","C")])
  K <- sum(singleNucFre[c("G","T")])
  M <- sum(singleNucFre[c("C","A")])
  Y <- sum(singleNucFre[c("C","T")])
  R <- sum(singleNucFre[c("G","A")])
  N <- 1
  nucFre <- c(singleNucFre,B,D,H,K,M,N,R,S,V,W,Y)
  names(nucFre) <- c(names(singleNucFre),"B","D","H","K","M","N","R","S","V","W","Y")
  ## get expected frequency of fasta
  expFreList <- do.call(append,lapply(seq.int(length(thispattern)),function(i){
    subPattern <- thispattern[i]
    rcSubPattern <- reverseComplement(DNAString(subPattern))
    patternTab <- table(unlist(strsplit(subPattern,"")))
    rcPatternTab <- table(unlist(strsplit(as.character(rcSubPattern),"")))
    forwardExpFre <- prod(nucFre[names(patternTab)]^patternTab)
    backwardExpFre <- prod(nucFre[names(rcPatternTab)]^rcPatternTab)
    names(forwardExpFre) <- names(backwardExpFre) <- names(thispattern[i])
    forwardExpFre+backwardExpFre
  }))
  ## get real frequency of input seq in peak region
  XNmatrix <- matrix(0,length(patternName),3)
  rownames(XNmatrix) <- patternName
  colnames(XNmatrix) <- c("x","n","prop.background")
  for(i in unique(width(oligoNset))){
    oligoNVec <- colSums(Biostrings::oligonucleotideFrequency(DNAStringSet(seqs$sequence),width=i))
    n <- which(width(oligoNset) %in% i)
    XNmatrix[n,1:2] <- do.call(rbind,lapply(n,function(j){
      pos.plus = gregexpr(thispattern.t[j], names(oligoNVec), perl = TRUE)
      pos.minus = gregexpr(revcomp.pattern.t[j], names(oligoNVec), perl = TRUE)
      c(sum(oligoNVec[which(sapply(pos.plus, "[[", 1) > 0)]) + sum(oligoNVec[which(sapply(pos.minus, "[[", 1) > 0)]),sum(oligoNVec))
    }))
  }
  XNmatrix[,3] <- expFreList
  binom.pvalue <- do.call(rbind,lapply(seq.int(nrow(XNmatrix)),function(i){
    binom.result <- stats::binom.test(XNmatrix[i,1],XNmatrix[i,2],XNmatrix[i,3],alternative="greater")
    binom.result$p.value
  }))
  XNmatrix <- data.frame(XNmatrix,binom.pvalue)
  return(XNmatrix)
}
