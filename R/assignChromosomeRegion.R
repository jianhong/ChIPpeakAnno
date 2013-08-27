#import GenomicFeatures
## Jianhong Ou @ Mar.20, 2013
assignChromosomeRegion <-
function(peaks.RD, exon, TSS, utr5, utr3, proximal.promoter.cutoff=1000L, immediate.downstream.cutoff=1000L, nucleotideLevel=FALSE, precedence=NULL, TranscriptDb=NULL)
{
  ##check inputs
  if(!is.null(TranscriptDb)){
    if(!inherits(TranscriptDb, "TranscriptDb")) 
      stop("TranscriptDb must be an object of TranscriptDb, try\n?TranscriptDb\tto see more info.")
    if(!inherits(peaks.RD, c("RangedData","GRanges"))) stop("peaks.RD must be a RangedData or GRanges object.")
    if(!is.null(precedence)) {
      if(!all(precedence %in% c("Exons", "Introns", "fiveUTRs", "threeUTRs", "Promoters", "immediateDownstream"))) stop("precedence must be a combination of Exons, Introns, fiveUTRs, threeUTRs, Promoters, immediateDownstream")
    }
    if(inherits(peaks.RD, "RangedData")) peaks.RD <- as(peaks.RD, "GRanges")
    exons <- exons(TranscriptDb, columns=NULL)
    introns <- unique(unlist(intronsByTranscript(TranscriptDb)))
    fiveUTRs <- unique(unlist(fiveUTRsByTranscript(TranscriptDb)))
    threeUTRs <- unique(unlist(threeUTRsByTranscript(TranscriptDb)))
    transcripts <- transcripts(TranscriptDb, columns=NULL)
    options(warn = -1)
    try({
      promoters <- promoters(TranscriptDb, upstream=proximal.promoter.cutoff, downstream=0)
      immediateDownstream <- flank(transcripts, width=immediate.downstream.cutoff, start=FALSE, use.names=FALSE)
    })
    microRNAs <- tryCatch(microRNAs(TranscriptDb), error=function(e) return(NULL))
    tRNAs <- tryCatch(tRNAs(TranscriptDb), error=function(e) return(NULL))
    options(warn = 0)
    annotation <- list(exons, introns, fiveUTRs, threeUTRs, promoters, immediateDownstream)
    if(!is.null(microRNAs)) annotation <- c(annotation, "microRNAs"=microRNAs)
    if(!is.null(tRNAs)) annotation <- c(annotation, "tRNAs"=tRNAs)
    annotation <- lapply(annotation, function(.anno){mcols(.anno)<-NULL; .anno})
    names(annotation)[1:6] <- c("Exons", "Introns", "fiveUTRs", "threeUTRs", "Promoters", "immediateDownstream")
    ###clear seqnames, the format should be chr+NUM
    formatSeqnames <- function(gr) {
      seqlevels(gr)[grepl("^(\\d+|MT|M|X|Y)$", seqlevels(gr))] <-
        paste("chr", seqlevels(gr)[grepl("^(\\d+|MT|M|X|Y)$", seqlevels(gr))], sep="")
      seqlevels(gr)[seqlevels(gr)=="chrMT"] <- "chrM"
      gr
    }
    peaks.RD <- formatSeqnames(peaks.RD)
    peaks.RD <- unique(peaks.RD)
    annotation <- lapply(annotation, formatSeqnames)
    annotation <- GRangesList(annotation)
    newAnno <- c(unlist(annotation))
    if(!all(seqlevels(peaks.RD) %in% seqlevels(newAnno))){
      warning("peaks.RD has sequence levels not in TranscriptDb.")
      sharedlevels <- intersect(seqlevels(newAnno), seqlevels(peaks.RD))
      peaks.RD <- keepSeqlevels(peaks.RD, sharedlevels)
    }
    mcols(peaks.RD) <- NULL
    if(!is.null(precedence)){
      annotation <- annotation[unique(c(precedence,names(annotation)))]
    }
    annotation$Intergenic.Region <- peaks.RD
    anno.names <- names(annotation)
    if(nucleotideLevel){
      ##create a new annotation GRanges
      newAnno <- c(newAnno, peaks.RD)
      ##create a splited cluster (no overlapps for all ranges)
      newAnno <- disjoin(newAnno)
      ol.anno <- findOverlaps(newAnno, annotation)
      ##calculate Jaccard index
      ol.anno.splited <- split(queryHits(ol.anno),anno.names[subjectHits(ol.anno)])
      Intergenic.Region <- ol.anno.splited$Intergenic.Region
      jaccardIndex <- unlist(lapply(ol.anno.splited, function(.ele){
        intersection <- intersect(.ele, Intergenic.Region)
        union <- union(.ele, Intergenic.Region)
        length(intersection)/length(union)
      }))
      jaccardIndex <- jaccardIndex[anno.names]
    }else{
      ol.anno <- findOverlaps(peaks.RD, annotation)
      ##calculate Jaccard index
      ol.anno.splited <- split(queryHits(ol.anno),anno.names[subjectHits(ol.anno)])
      Intergenic.Region <- ol.anno.splited$Intergenic.Region
      jaccardIndex <- unlist(lapply(anno.names, function(.name){
        union <- length(Intergenic.Region)+length(annotation[[.name]])-length(ol.anno.splited[[.name]])
        intersection <- length(ol.anno.splited[[.name]])
        intersection/union
      }))
      names(jaccardIndex) <- anno.names
    }
    ol.anno <- as.data.frame(ol.anno)
    if(nucleotideLevel){
      ####keep the part in peaks.RD
      ol.anno <- ol.anno[ol.anno[,1] %in% unique(ol.anno[ol.anno[,2]==length(annotation),1]),]
    }
    ####keep the part only annotated in peaks.RD for Intergenic.Region
    ol.anno.splited <- split(ol.anno, ol.anno[,2])
    hasAnnoHits <- do.call(rbind, ol.anno.splited[names(ol.anno.splited)!=as.character(length(anno.names))])
    hasAnnoHits <- unique(hasAnnoHits[,1])
    ol.anno <- ol.anno[!(ol.anno[,2]==length(anno.names) & (ol.anno[,1] %in% hasAnnoHits)), ]    
    ###recalculate Jaccard index for Intergenic.Region
    jaccardIndex["Intergenic.Region"] <- nrow(ol.anno[ol.anno[,2]==as.character(length(anno.names)),])/length(Intergenic.Region)
    if(!is.null(precedence)){
      ol.anno <- ol.anno[!duplicated(ol.anno[,1]), ]
    }
    ##calculate percentage
    if(nucleotideLevel){
      queryHits <- newAnno[ol.anno[,1]]
      totalLen <- sum(width(queryHits[!duplicated(queryHits)]))
      ##mcols(queryHits)$subjectHits <- anno.names[ol.anno[,2]]
      queryHits.list <- split(queryHits, anno.names[ol.anno[,2]])
      lens <- unlist(lapply(queryHits.list, function(.ele) sum(width(.ele))))
      percentage <- 100 * lens/totalLen
    }else{
      subjectHits <-anno.names[ol.anno[,2]]
      counts <- table(subjectHits)
      percentage <- 100 * counts / length(peaks.RD)
    }
    len <- length(anno.names) - length(percentage)
    if(len>0) {
      tobeadd <- rep(0, len)
      names(tobeadd) <- anno.names[!anno.names %in% names(percentage)]
      percentage <- c(percentage, tobeadd)
    }
    percentage <- percentage[anno.names]
    return(list(percentage=percentage, jaccard=jaccardIndex))
  }else{
    message("Please try to use TranscriptDb next time. Try\n?TranscriptDb\tto see more info.")
    annotationList <- list(exon, TSS, utr5, utr3)
    names(annotationList) <- c("Exon", "TSS", "UTR5", "UTR3")
    status <- lapply(annotationList, function(.ele) {
      if(!inherits(.ele, "RangedData")){
        stop("Annotation of exon, TSS, utr5, utr3 must be objects of RangedData.")
      }
    })
    if(!inherits(peaks.RD, "RangedData")) stop("peaks.RD must be a RangedData object.") 
    ann.peaks = annotatePeakInBatch(peaks.RD, AnnotationData = TSS)
    ann.peaks = ann.peaks[!is.na(ann.peaks$distancetoFeature),]
    upstream = ann.peaks[ ann.peaks$insideFeature=="upstream" | (ann.peaks$distancetoFeature<0 & ann.peaks$insideFeature =="overlapStart" & abs(ann.peaks$distancetoFeature)>ann.peaks$shortestDistance ) | ann.peaks$insideFeature=="includeFeature" | (ann.peaks$distancetoFeature>=0 & ann.peaks$insideFeature =="overlapStart" & ann.peaks$distancetoFeature ==ann.peaks$shortestDistance),]
    
    proximal.promoter.n = dim(upstream[upstream$distancetoFeature>= - proximal.promoter.cutoff | upstream$shortestDistance<= proximal.promoter.cutoff, ])[1]
    enhancer.n = dim(upstream)[1] - proximal.promoter.n
    
    downstream = ann.peaks[ann.peaks$insideFeature =="downstream",]
    immediateDownstream.n = dim(downstream[downstream$distancetoFeature <= immediate.downstream.cutoff,])[1]
    enhancer.n = enhancer.n +  dim(downstream[downstream$distancetoFeature > immediate.downstream.cutoff,])[1]
    
    inside.peaks =  ann.peaks[ann.peaks$insideFeature =="inside" | ann.peaks$insideFeature =="overlapEnd" |  (ann.peaks$insideFeature =="overlapStart" & ann.peaks$distancetoFeature >=0 & ann.peaks$distancetoFeature != ann.peaks$shortestDistance) | (ann.peaks$insideFeature =="overlapStart" & ann.peaks$distancetoFeature <0 & abs(ann.peaks$distancetoFeature) ==ann.peaks$shortestDistance) ,]
    
    ann.utr5.peaks = annotatePeakInBatch(inside.peaks, AnnotationData = utr5)
    
    proximal.promoter.n = proximal.promoter.n + dim(ann.utr5.peaks[ann.utr5.peaks$insideFeature =="upstream",])[1]
    
    utr5.n = dim(ann.utr5.peaks[ann.utr5.peaks$insideFeature %in% c("includeFeature" , "inside") | (ann.utr5.peaks$insideFeature =="overlapStart"  & ann.utr5.peaks$distancetoFeature >=0 & ann.utr5.peaks$distancetoFeature != ann.utr5.peaks$shortestDistance)  | (ann.utr5.peaks$insideFeature =="overlapStart" & ann.utr5.peaks$distancetoFeature <0 & abs(ann.utr5.peaks$distancetoFeature) ==ann.utr5.peaks$shortestDistance)  | (ann.utr5.peaks$insideFeature =="overlapEnd" & ann.utr5.peaks$strand=="+" & abs(start(ann.utr5.peaks)-ann.utr5.peaks$end_position) >= (end(ann.utr5.peaks)-ann.utr5.peaks$end_position)) | (ann.utr5.peaks$insideFeature =="overlapEnd" & ann.utr5.peaks$strand=="-" & abs(end(ann.utr5.peaks)-ann.utr5.peaks$start_position) >= abs(start(ann.utr5.peaks)-ann.utr5.peaks$start_position )), ] )[1]
    
    proximal.promoter.n = proximal.promoter.n +  dim(ann.utr5.peaks[(ann.utr5.peaks$insideFeature =="overlapStart"  & ann.utr5.peaks$distancetoFeature >=0 & ann.utr5.peaks$distancetoFeature == ann.utr5.peaks$shortestDistance)  | (ann.utr5.peaks$insideFeature =="overlapStart" & ann.utr5.peaks$distancetoFeature <0 & abs(ann.utr5.peaks$distancetoFeature) !=ann.utr5.peaks$shortestDistance),])[1]
    
    downstream.utr5 = ann.utr5.peaks[ann.utr5.peaks$insideFeature =="downstream" | (ann.utr5.peaks$insideFeature =="overlapEnd" & ann.utr5.peaks$strand=="+" & abs(start(ann.utr5.peaks)-ann.utr5.peaks$end_position) < (end(ann.utr5.peaks)-ann.utr5.peaks$end_position)) | (ann.utr5.peaks$insideFeature =="overlapEnd" & ann.utr5.peaks$strand=="-" & abs(end(ann.utr5.peaks)-ann.utr5.peaks$start_position) < abs(start(ann.utr5.peaks)-ann.utr5.peaks$start_position )), ] 
    
    ann.utr3.peaks = annotatePeakInBatch(downstream.utr5, AnnotationData = utr3)
    
    utr3.n = dim(ann.utr3.peaks[ann.utr3.peaks$insideFeature %in% c("includeFeature" , "overlapStart", "overlapEnd", "inside"),])[1]
    
    rest.peaks = ann.utr3.peaks[ann.utr3.peaks$insideFeature %in% c("downstream", "upstream"),]
    
    ann.rest.peaks = annotatePeakInBatch(rest.peaks, AnnotationData = exon)
    
    intron.n = dim(ann.rest.peaks[ann.rest.peaks$insideFeature %in% c("downstream", "upstream"),])[1]
    exon.n =dim(ann.rest.peaks)[1] - intron.n
    
    total = dim(peaks.RD)[1]/100
    
    list( "Exons" =exon.n/total, "Introns"=intron.n/total, "fiveUTRs" = utr5.n/total, "threeUTRs" = utr3.n/total, "Promoters"= proximal.promoter.n/total, "immediate.Downstream" = immediateDownstream.n/total, "Intergenic.Region" = enhancer.n/total)
  }
}
