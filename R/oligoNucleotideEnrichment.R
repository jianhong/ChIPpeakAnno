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
#' @param upstream upstream of peak
#' 
#' @param downstream downstream of peak
#'
#' @param genome BSgenome object or mart object. Please refer to available.genomes in BSgenome package and useMart in bioMaRt package for details
#'
#' @param methodBackground The method to get the background of compared oligonucleotide. "chromSelectRandomly" (default) is used to select background chromosomes from all chromosomesor "shuffle".
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
#' Please see \link[universalmotif]{shuffle_sequences} for the paramter k, method, nthreads, window, window.size and window.overlap.
#'
#' @author Junhui Li
#'
#' @importFrom  Biostrings oligonucleotideFrequency
#' 
#' @importFrom  Biostrings readDNAStringSet
#' 
#' @importFrom  Biostrings reverseComplement
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
#' background="chromSelectRandomly")
#' }
#' @export

oligoNucleotideEnrichment <- function(filepath,
                                      format="fasta",
                                      peaks,
                                      upstream = 0,
                                      downstream = 0,
                                      genome,
                                      background = c("chromSelectRandomly","shuffle"),
                                      chromosome = NULL,
                                      ...,
                                      times = 1000,
                                      alpha = 0.05){
  stopifnot("file doesn't exist!"=file.exists(filepath))
  background <- match.arg(methodBackground)
  stopifnot("genome is required parameter, 
           please pass in either a BSgenome object or a Mart object."=
              is(genome, "BSgenome") | is(genome, "Mart"))
  stopifnot("The 'times' parameter should be increased to a sufficient extent."=round(times*alpha,0) != 0)
  nPeak <- length(peaks)
  stopifnot("There is no peak, please check your data."=nPeak != 0)
  allseqLen <- seqlengths(genome)
  seqWidth <- width(peaks)
  candiChrom <- names(allseqLen >= max(seqWidth))
  stopifnot("The length of peak sequence should be less than the length of chromosomes"=
              all(seqWidth < allseqLen[as.character(seqnames(peaks)@values)]))
  peakSeq <- getAllPeakSequence(peaks, upstream = upstream, downstream = downstream, genome = genome)
  peakPvalue <- oligoNucleotideSummary(filepath, format = format, seqs = peakSeq)
  
  backgroundStatDistribution <- matrix(NA,0,4)
  colnames(backgroundStatDistribution) <- c("x","n","prop_background","binom_pvalue")
  if(methodBackground == "chromSelectRandomly"){
    backgroundStatDistribution <- do.call(rbind,lapply(seq.int(times),function(n) {
      if(is.null(chromosome)){
        chrs <- as.character(seqnames(peaks))
        givenSeqLen <- allseqLen[chrs]
      }else{
        chrs <- sample(candiChrom,nPeak)
        givenSeqLen <- allseqLen[chrs]
      }
      startPos <- vapply(givenSeqLen-seqWidth,function(x){sample(seq.int(x),1)},numeric(1))
      endPos <- startPos + seqWidth - 1
      backGroundPeak <- GRanges(seqnames = chrs,
                          IRanges(start = startPos,
                                  end = endPos,
                                  names = paste0("peak",seq.int(nPeak))))
      backGroundSeq <- getAllPeakSequence(backGroundPeak, upstream = 0, downstream = 0, genome = genome)
      statEachTime <- oligoNucleotideSummary(filepath,format=format,seqs=backGroundSeq)
      rownames(statEachTime) <- paste0(n,"_",rownames(statEachTime))
      statEachTime
    }))
  }else if(methodBackground=="shuffle"){
    colnames(mcols(peakSeq))[3] <- "orig.sequence"
    names(peakSeq$orig.sequence) <- do.call(paste, c(as.data.frame(peakSeq)[,1:3], sep="_"))
    inputSeq <- DNAStringSet(peakSeq$orig.sequence)
    backgroundStatDistribution <- do.call(rbind,lapply(seq.int(times),function(n){
      sequences.shuffled <- universalmotif::shuffle_sequences(inputSeq, ...)
      peakSeq$sequence <- as.vector(sequences.shuffled)
      statEachTime <- oligoNucleotideSummary(filepath,format=format,seqs=peakSeq)
      rownames(statEachTime) <- paste0(n,"_",rownames(statEachTime))
      statEachTime
    }))
  }

  backgroundPvalue <- do.call(cbind,lapply(seq.int(nrow(peakPvalue)),function(i){
    sort(as.numeric(backgroundStatDistribution[seq(i,nrow(backgroundStatDistribution),nrow(peakPvalue)),4]))
  }))
  colnames(backgroundPvalue) <- seq.int(nrow(peakPvalue))
  peakPvalue$cutoff <- backgroundPvalue[round(times*alpha,0),]
  return(peakPvalue)
}

oligoNucleotideSummary <- function(filepath,format="fasta",seqs){
  stopifnot("File doesn't exist!"=file.exists(filepath))
  oligoNset = readDNAStringSet(filepath=filepath, format=format, use.names = TRUE)
  patternName <- names(oligoNset)
  revcomp.pattern = reverseComplement(oligoNset)
  thispattern = as.character(oligoNset)
  thispattern.t = translatePattern(thispattern)
  revcomp.pattern = as.character(revcomp.pattern)
  revcomp.pattern.t = translatePattern(revcomp.pattern)

  ## get frequency of single nuleotide in seqs
  ACGTcount <- colSums(oligonucleotideFrequency(DNAStringSet(seqs$sequence),width=1))
  allBaseFreq <- ACGTcount/sum(ACGTcount)
  allBaseFreq['D'] <- sum(ACGTfreq[c("G","A","T")])
  allBaseFreq['H'] <- sum(ACGTfreq[c("C","A","T")])
  allBaseFreq['B'] <- sum(ACGTfreq[c("G","C","T")])
  allBaseFreq['V'] <- sum(ACGTfreq[c("G","A","C")])
  allBaseFreq['W'] <- sum(ACGTfreq[c("A","T")])
  allBaseFreq['S'] <- sum(ACGTfreq[c("G","C")])
  allBaseFreq['K'] <- sum(ACGTfreq[c("G","T")])
  allBaseFreq['M'] <- sum(ACGTfreq[c("C","A")])
  allBaseFreq['Y'] <- sum(ACGTfreq[c("C","T")])
  allBaseFreq['R'] <- sum(ACGTfreq[c("G","A")])
  allBaseFreq['N'] <- 1
  
  ## get expected frequency of fasta
  expFreqSet <- sapply(thispattern,function(i){
    rcSubPattern <- toString(reverseComplement(DNAString(i)))
    patternTab <- table(strsplit(i,"")[[1]])
    rcPatternTab <- table(strsplit(rcSubPattern,"")[[1]])
    forwardExpFre <- prod(allBaseFreq[names(patternTab)]^patternTab)
    backwardExpFre <- prod(allBaseFreq[names(rcPatternTab)]^rcPatternTab)
    expFreq <- forwardExpFre+backwardExpFre
    names(expFreq) <- names(i)
    expFreq
  })

  ## get real frequency of input seq in peak region
  peakFreq <- matrix(0,length(patternName),3)
  rownames(peakFreq) <- patternName
  colnames(peakFreq) <- c("x","n","prop_background")
  peakFreq[,3] <- expFreqSet
  
  peakFreq[,1:2] <- t(sapply(seq.int(width(oligoNset)),function(i){
    oligoNVec <- colSums(oligonucleotideFrequency(DNAStringSet(seqs$sequence),width=width(oligoNset)[i]))
    pos.plus = gregexpr(thispattern.t[i], names(oligoNVec), perl = TRUE)
    pos.minus = gregexpr(revcomp.pattern.t[i], names(oligoNVec), perl = TRUE)
    c(sum(oligoNVec[which(sapply(pos.plus, "[[", 1) > 0)]) + sum(oligoNVec[which(sapply(pos.minus, "[[", 1) > 0)]),sum(oligoNVec))
  }))

  binom_pvalue <- sapply(seq.int(nrow(peakFreq)),function(i){
    stats::binom.test(peakFreq[i,1],peakFreq[i,2],peakFreq[i,3],alternative="greater")$p.value
  })
  peakFreq <- data.frame(peakFreq,binom_pvalue)
  return(peakFreq)
}
