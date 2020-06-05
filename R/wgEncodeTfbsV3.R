#' transcription factor binding site clusters (V3) from ENCODE
#' 
#' possible binding pool for human (hg19) from transcription factor binding
#' site clusters (V3) from ENCODE data and removed the HOT spots
#' 
#' How to generate the data:
#' 
#' temp <- tempfile()
#' 
#' download.file(file.path("http://hgdownload.cse.ucsc.edu", "goldenPath",
#' 
#' "hg19", "encodeDCC",
#' 
#' "wgEncodeRegTfbsClustered",
#' 
#' "wgEncodeRegTfbsClusteredV3.bed.gz"), temp)
#' 
#' data <- read.delim(gzfile(temp, "r"), header=FALSE)
#' 
#' unlink(temp)
#' 
#' colnames(data)[1:4] <- c("seqnames", "start", "end", "TF")
#' 
#' wgEncodeRegTfbsClusteredV3 <- GRanges(as.character(data$seqnames),
#' 
#' IRanges(data$start, data$end),
#' 
#' TF=data$TF)
#' 
#' data(HOT.spots)
#' 
#' hot <- reduce(unlist(HOT.spots))
#' 
#' ol <- findOverlaps(wgEncodeRegTfbsClusteredV3, hot)
#' 
#' wgEncodeTfbsV3 <- wgEncodeRegTfbsClusteredV3[-unique(queryHits(ol))]
#' 
#' wgEncodeTfbsV3 <- reduce(wgEncodeTfbsV3)
#' 
#' save(list="wgEncodeTfbsV3",
#' 
#' file="data/wgEncodeTfbsV3.rda",
#' 
#' compress="xz", compression_level=9)
#' 
#' @name wgEncodeTfbsV3
#' @docType data
#' @format An object of GRanges.
#' @source http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/
#' wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV3.bed.gz
#' @keywords datasets
#' @examples
#' 
#' data(wgEncodeTfbsV3)
#' head(wgEncodeTfbsV3)
#' 
"wgEncodeTfbsV3"
