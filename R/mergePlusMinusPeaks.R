#' Merge peaks from plus strand and minus strand
#' 
#' Merge peaks from plus strand and minus strand within certain distance apart,
#' and output merged peaks as bed format.
#' 
#' 
#' @param peaks.file Specify the peak file. The peak file should contain peaks
#' from both plus and minus strand
#' @param columns Specify the column names in the peak file
#' @param sep Specify column delimiter, default tab-delimited
#' @param header Specify whether the file has a header row, default TRUE
#' @param distance.threshold Specify the maximum gap allowed between the plus
#' stranded and the nagative stranded peak
#' @param plus.strand.start.gt.minus.strand.end Specify whether plus strand
#' peak start greater than the paired negative strand peak end. Default to TRUE
#' @param output.bedfile Specify the bed output file name
#' @return output the merged peaks in bed file and a data frame of the bed
#' format
#' @author Lihua Julie Zhu
#' @seealso annotatePeakInBatch, findOverlappingPeaks, makeVennDiagram
#' @references Zhu L.J. et al. (2010) ChIPpeakAnno: a Bioconductor package to
#' annotate ChIP-seq and ChIP-chip data. BMC Bioinformatics 2010,
#' 11:237doi:10.1186/1471-2105-11-237
#' @keywords misc
#' @export
#' @importFrom matrixStats rowMins rowMaxs
#' @importFrom utils read.table write.table
#' @examples
#' 
#' 
#' if (interactive())
#' {
#'     data(myPeakList)
#'     data(TSS.human.NCBI36)
#'     library(matrixStats)
#'         peaks <- system.file("extdata", "guide-seq-peaks.txt", 
#'                               package = "ChIPpeakAnno")
#'         merged.bed <- mergePlusMinusPeaks(peaks.file = peaks, 
#'                                           columns=c("name", "chromosome", 
#'                                                     "start", "end", "strand", 
#'                                                     "count", "count"), 
#'                                           sep = "\t", header = TRUE,  
#'                                           distance.threshold = 100,  
#'                                 plus.strand.start.gt.minus.strand.end = TRUE, 
#'                                           output.bedfile = "T2test100bp.bed")
#' }
#' 
mergePlusMinusPeaks <- function(peaks.file, 
                                columns=c("name", "chromosome", "start", 
                                          "end", "strand", "count", 
                                          "count", "count", "count"), 
                                sep = "\t", header = TRUE,  
                                distance.threshold = 100,  
                                plus.strand.start.gt.minus.strand.end = TRUE, 
                                output.bedfile)
{
    peaks<- read.table (peaks.file, sep = sep, header = header)
    if (dim(peaks)[2] != length(columns))
        stop("Number of columns specified differs from the number of columns 
             in the input peak file, please modify the column specification 
             in parameter columns accordingly!")
    if (length(intersect(columns, "strand")) != 1)
        stop("Please include one strand column in columns and 
             corresponding peaks.file")
    colnames(peaks) <- columns
    pos.peaks <- subset(peaks, peaks[, which(columns == "strand")] == "+")
    neg.peaks <- subset(peaks, peaks[, which(columns == "strand")] == "-")
    if (dim(pos.peaks)[1] == 0 || dim(neg.peaks)[1] == 0)
        stop("Need peaks from both + and - strand to merge")
    pos.gr <- makeGRangesFromDataFrame(pos.peaks)
    neg.gr <- makeGRangesFromDataFrame(neg.peaks)
    names(pos.gr) <- paste(paste(seqnames(pos.gr), strand(pos.gr), sep=""),
                           start(pos.gr),end(pos.gr), sep=":")
    names(neg.gr) <- paste(paste(seqnames(neg.gr), strand(neg.gr), sep=""),
                           start(neg.gr),end(neg.gr), sep=":")
    if (plus.strand.start.gt.minus.strand.end)
        gr <- annotatePeakInBatch(pos.gr, featureType = "TSS",
                                  AnnotationData = neg.gr, 
                                  output="nearestStart", 
                                  PeakLocForDistance = "TSS")
    else
        gr <- annotatePeakInBatch(neg.gr, featureType = "TSS", 
                                  AnnotationData = pos.gr, 
                                  output="nearestStart",
                                  PeakLocForDistance = "end")
    pos.peaks <- 
        cbind(peak = paste(paste(pos.peaks[,which(columns == "chromosome")], 
                                 "+", sep = ""), 
                           pos.peaks[,which(columns == "start")], 
                           pos.peaks[,which(columns == "end")], sep=":"), 
              pos.peaks)
    neg.peaks <- 
        cbind(feature = paste(paste(neg.peaks[,which(columns == "chromosome")],
                                    "-", sep = ""), 
                              neg.peaks[,which(columns == "start")], 
                              neg.peaks[,which(columns == "end")], sep=":"), 
              neg.peaks)
    if (!plus.strand.start.gt.minus.strand.end)
    {
        colnames(pos.peaks)[1] <-  "feature"
        colnames(neg.peaks)[1] <-  "peak"
    }
    new.count.columns <- which(columns == "count") + 1
    colnames(neg.peaks)[new.count.columns] <- 
        paste("minus", colnames(neg.peaks)[new.count.columns], sep=":")
    colnames(pos.peaks)[new.count.columns] <- 
        paste("plus", colnames(pos.peaks)[new.count.columns], sep=":")
    pos.peaks <- pos.peaks[, c(1, new.count.columns)]
    neg.peaks <- neg.peaks[, c(1, new.count.columns)]
    ann.peaks <- as.data.frame(gr[abs(gr$distancetoFeature) <= 
                                      distance.threshold & 
                                      gr$distancetoFeature < 0,])
    temp <- merge(pos.peaks, ann.peaks)
    temp1 <- merge(neg.peaks, temp)
    temp1 <- cbind(temp1, 
                   totalCount = rowSums(temp1[,grep("count", 
                                                    colnames(temp1))]))
    temp1$names = paste(temp1$peak, temp1$feature, sep=":")
    temp1$minStart <- rowMins(as.matrix(temp1[, c("start_position", "start")]))
    temp1$maxEnd <- rowMaxs(as.matrix(temp1[, c("end_position", "end")]))
    bed.temp <- 
        temp1[, c("seqnames", "minStart", "maxEnd", "names", "totalCount")]
    bed.temp <- cbind(bed.temp, strand = "+")
    write.table(bed.temp, file = output.bedfile, sep = "\t", col.names = FALSE,
                row.names=FALSE, quote=FALSE)
    bed.temp
}
