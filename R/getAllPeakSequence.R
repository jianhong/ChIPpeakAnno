#' Obtain genomic sequences around the peaks
#' 
#' Obtain genomic sequences around the peaks leveraging the BSgenome and
#' biomaRt package
#' 
#' 
#' @param myPeakList An object of \link[GenomicRanges:GRanges-class]{GRanges}:
#' See example below
#' @param upstream upstream offset from the peak start, e.g., 200
#' @param downstream downstream offset from the peak end, e.g., 200
#' @param genome BSgenome object or mart object. Please refer to
#' available.genomes in BSgenome package and useMart in bioMaRt package for
#' details
#' @param AnnotationData GRanges object with annotation information.
#' @return \link[GenomicRanges:GRanges-class]{GRanges} with slot start holding
#' the start position of the peak, slot end holding the end position of the
#' peak, slot rownames holding the id of the peak and slot seqnames holding the
#' chromosome where the peak is located. In addition, the following variables
#' are included: \item{upstream}{upstream offset from the peak start}
#' \item{downstream}{downstream offset from the peak end} \item{sequence}{the
#' sequence obtained}
#' @author Lihua Julie Zhu, Jianhong Ou
#' @references Durinck S. et al. (2005) BioMart and Bioconductor: a powerful
#' link between biological biomarts and microarray data analysis.
#' Bioinformatics, 21, 3439-3440.
#' @keywords misc
#' @export
#' @importFrom BiocGenerics start end width strand
#' @importFrom GenomeInfoDb seqlengths seqlevels seqnames `seqlengths<-`
#' @importFrom Biostrings getSeq
#' @examples
#' 
#' #### use Annotation data from BSgenome
#' peaks <- GRanges(seqnames=c("NC_008253", "NC_010468"),
#'                  IRanges(start=c(100, 500), end=c(300, 600), 
#'                          names=c("peak1", "peak2")))
#' library(BSgenome.Ecoli.NCBI.20080805)
#' seq <- getAllPeakSequence(peaks, upstream=20, downstream=20, genome=Ecoli)
#' write2FASTA(seq, file="test.fa")
#' 
getAllPeakSequence <- function(myPeakList, 
                               upstream=200L, downstream=upstream, 
                               genome, AnnotationData)
{
    if (!inherits(myPeakList, c("GRanges"))) {
        stop("No valid myPeakList passed in. It needs to be GRanges object")
    }
    if (missing(genome))
    {
        stop("genome is required parameter, 
             please pass in either a BSgenome object or a Mart object!")
    }
  old_name <- names(myPeakList)
    if(length(old_name)!=length(myPeakList)){
      names(myPeakList) <- paste0("peak_", seq_along(myPeakList))
    }
    myPeakList.bk <- myPeakList
    if (is(genome, "BSgenome"))
    {
        strand <- strand(myPeakList)
        width <- width(myPeakList)
        ##  start(myPeakList)
        start(myPeakList) <- ifelse(strand=="-", 
                                    start(myPeakList) - as.numeric(downstream),
                                    start(myPeakList) - as.numeric(upstream))
        lt0 <- start(myPeakList)<=0
        start <- start(myPeakList)-1
        start[!lt0] <- 0
        start(myPeakList)[lt0] <- 1
        width(myPeakList) <- width + as.numeric(upstream) + 
            as.numeric(downstream) + start
        strand[strand!="-"] <- "+"
        strand(myPeakList) <- strand
        ##how to make the seqname be same is a question.
        if(!all(seqlevels(myPeakList) %in% seqnames(genome))){
            myPeakList <- formatSeqnames(myPeakList)
            genome <- formatSeqnames(genome)
        }
        ends <- unlist(
            apply(
                cbind(end(myPeakList), 
                      seqlengths(genome)[as.character(seqnames(myPeakList))]), 
                1,
                min))
        ends <- unname(ifelse(is.na(ends), end(myPeakList), ends))
        keep <- start(myPeakList) <= ends
        end(myPeakList)[keep] <- ends[keep]
        myPeakList <- myPeakList[keep]
        seq <- getSeq(genome, myPeakList, as.character=TRUE)
        
        myPeakList <- myPeakList.bk
        myPeakList$upstream <- rep(upstream,length(myPeakList.bk))
        myPeakList$downstream <- rep(downstream,length(myPeakList.bk))
        myPeakList$sequence <- NA
        if(!all(names(myPeakList) %in% names(seq))){
          warning("The genome assembly are not identical to your peaks!")
        }
        myPeakList[names(seq)]$sequence <- seq
        if(all(seqlevels(myPeakList) %in% seqlevels(genome))){
            seqlengths(myPeakList) <- seqlengths(genome)[seqlevels(myPeakList)]
        }   
        names(myPeakList) <- old_name
        myPeakList
    }else if (is(genome, "Mart")){
        if (missing(AnnotationData)) {
            message("No AnnotationData as GRanges is passed in, 
                    so now querying biomart database for AnnotationData ....")
            AnnotationData <- getAnnotation(genome)
            message("Done querying biomart database, start annotating ....
                    Better way would be calling getAnnotation before 
                    querying for sequence")
        }
        if (!inherits(AnnotationData, c("GRanges"))) {
            stop("AnnotationData needs to be GRanges object. 
                 Better way would be calling getAnnotation.")
        }
        
        downstream.bk = downstream
        plusAnno = AnnotationData[strand(AnnotationData)=="+"]
        temp = annotatePeakInBatch(myPeakList, AnnotationData=plusAnno)
        TSSlength =temp$end_position - temp$start_position
        downstream = end(temp) - start(temp) + downstream
        temp$downstream=  downstream
        temp$TSSlength = TSSlength
        
        myList3 = as.data.frame(temp)
        rm(temp)
        
        if (dim(myList3)[1] != 0)
        {
            l3 =  cbind(as.character(myList3$feature), 
                        as.numeric(as.character(myList3$distancetoFeature)), 
                        as.numeric(rep(upstream, dim(myList3)[1])), 
                        as.numeric(as.character(myList3$downstream)), 
                        as.numeric(as.character(myList3$start_position)), 
                        as.numeric(as.character(myList3$end_position)))
            r3 = apply(l3, 1, getGeneSeq, genome)
        }
        else
        {
            r3 = 0
        }
        
        if (is.list(r3))
        {
            r = as.data.frame(do.call("rbind",r3))
        }
        else
        {
            stop("No sequence found error!")
        }
        colnames(r)= c("feature", "distancetoFeature", "upstream", 
                       "downstream", "seq")
        r4 = merge(r, myList3)
        GRanges(seqnames=as.character(r4$space),
                ranges=IRanges(start=r4$start, 
                               end  =r4$end, 
                               names=as.character(r4$name)),
                strand="+",
                upstream=rep(upstream, dim(r4)[1]), 
                downstream= rep(downstream.bk, dim(r4)[1]),
                sequence=unlist(r4$seq))
    }else{
        stop("genome needs to be either a BSgenome object or Mart object!")
    }
}
