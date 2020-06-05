#' Find possible enhancers depend on DNA interaction data
#' 
#' Find possible enhancers by data from chromosome conformation capture
#' techniques such as 3C, 5C or HiC.
#' 
#' 
#' @param peaks peak list, \link[GenomicRanges:GRanges-class]{GRanges} object
#' @param annoData annotation data, \link[GenomicRanges:GRanges-class]{GRanges}
#' object
#' @param DNAinteractiveData DNA interaction data,
#' \link[GenomicRanges:GRanges-class]{GRanges} object with interaction blocks
#' informations.
#' @param bindingType Specifying the criteria to associate peaks with
#' annotation. Here is how to use it together with the parameter bindingRegion.
#' The annotation will be shift to a new position depend on the DNA interaction
#' region.  \itemize{ \item To obtain peaks within 5kb upstream and up to 3kb
#' downstream of shift TSS within the gene body, set bindingType = "startSite"
#' and bindingRegion = c(-5000, 3000) \item To obtain peaks up to 5kb upstream
#' within the gene body and 3kb downstream of shift gene/Exon End, set
#' bindingType = "endSite" and bindingRegion = c(-5000, 3000) \item To obtain
#' peaks with nearest bi-directional enhancer regions within 5kb upstream and
#' 3kb downstream of shift TSS, set bindingType =
#' "nearestBiDirectionalPromoters" and bindingRegion = c(-5000, 3000) }
#' \describe{ \item{startSite}{start position of the feature (strand is
#' considered)} \item{endSite}{end position of the feature (strand is
#' considered)} \item{nearestBiDirectionalPromoters}{nearest enhancer regions
#' from both direction of the peaks (strand is considered). It will report
#' bidirectional enhancer regions if there are enhancer regions in both
#' directions in the given region (defined by bindingRegion). Otherwise, it
#' will report the closest enhancer regions in one direction.} }
#' @param bindingRegion Annotation range used together with bindingType, which
#' is a vector with two integer values, default to c (-5000, 5000). The first
#' one must be no bigger than 0. And the sec ond one must be no less than 1.
#' For details, see bindingType.
#' @param ignore.peak.strand ignore the peaks strand or not.
#' @param ...  Not used.
#' @return Output is a GRanges object of the annotated peaks.
#' @author Jianhong Ou
#' @seealso See Also as \code{\link{annotatePeakInBatch}}
#' @keywords misc
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom BiocGenerics start end width strand
#' @importFrom S4Vectors elementNROWS mcols queryHits subjectHits
#' @examples
#' 
#'   bed <- system.file("extdata", 
#'                      "wgEncodeUmassDekker5CGm12878PkV2.bed.gz",
#'                      package="ChIPpeakAnno")
#'   DNAinteractiveData <- toGRanges(gzfile(bed))
#'   library(EnsDb.Hsapiens.v75)
#'   annoData <- toGRanges(EnsDb.Hsapiens.v75, feature="gene")
#'   data("myPeakList")
#'   findEnhancers(myPeakList[500:1000], annoData, DNAinteractiveData)
#' 
findEnhancers <- function(peaks, annoData, DNAinteractiveData, 
                    bindingType=c("nearestBiDirectionalPromoters",
                                  "startSite", "endSite"), 
                    bindingRegion=c(-5000, 5000), 
                    ignore.peak.strand=TRUE, ...){
    stopifnot(inherits(peaks, "GRanges"))
    stopifnot(inherits(annoData, c("annoGR", "GRanges")))
    stopifnot(length(intersect(seqlevelsStyle(peaks),
                               seqlevelsStyle(annoData)))>0)
    stopifnot(inherits(DNAinteractiveData, "GRanges"))
    stopifnot(length(DNAinteractiveData$blocks)>0)
    stopifnot(all(elementNROWS(DNAinteractiveData$blocks)==2))
    stopifnot(length(intersect(seqlevelsStyle(peaks),
                               seqlevelsStyle(DNAinteractiveData)))>0)
    bindingType <- match.arg(bindingType)
    stopifnot(length(bindingRegion)==2)
    stopifnot(bindingRegion[1]<=0 && bindingRegion[2]>=1)
    
    if(inherits(annoData, "annoGR")){
        annoData <- as(annoData, "GRanges")
    }
    
    peaks$peak.oid.to.be.deleted <- 1:length(peaks)
    
    HiC_FOR <- lapply(DNAinteractiveData$blocks, `[`, 1)
    HiC_REV <- lapply(DNAinteractiveData$blocks, `[`, 2)
    HiC_FOR <- unlist(IRangesList(HiC_FOR))
    HiC_REV <- unlist(IRangesList(HiC_REV))
    ## BED file blocks are half open half close.
    HiC_FOR <- shift(HiC_FOR, start(DNAinteractiveData)-1)
    HiC_REV <- shift(HiC_REV, start(DNAinteractiveData)-1)
    HiC_FOR_GR <- HiC_REV_GR <- DNAinteractiveData
    ranges(HiC_FOR_GR) <- HiC_FOR
    ranges(HiC_REV_GR) <- HiC_REV
    
    HiC.A.ups <- shift(HiC_FOR_GR, -max(abs(bindingRegion)))
    width(HiC.A.ups) <- max(abs(bindingRegion))
    HiC.D.dws <- shift(HiC_REV_GR, max(abs(bindingRegion)))
    start(HiC.D.dws) <- end(HiC_REV_GR) + 1
    HiC.BC <- HiC_REV_GR
    start(HiC.BC) <- end(HiC_FOR_GR) + 1
    end(HiC.BC) <- start(HiC_REV_GR) - 1
    HiC.groups <- list(A.ups=HiC.A.ups,
                       AB=HiC_FOR_GR,
                       BC=HiC.BC,
                       CD=HiC_REV_GR,
                       D.dws=HiC.D.dws)
    peaks.ol.HiCdata <- lapply(HiC.groups, function(hic){
        ol <- findOverlaps(peaks, hic)
        this.peaks <- peaks[queryHits(ol)]
        this.peaks$HiC.idx <- subjectHits(ol)
        this.peaks
    })
    if(all(elementNROWS(peaks.ol.HiCdata)==0)){
        return(GRanges())
    }
    ## refine annoData by HiCdata
    anno.ol.HiCdata <- lapply(HiC.groups, function(hic){
        ol <- findOverlaps(annoData, hic)
        annoData.ol.HiC <- annoData[queryHits(ol)]
        annoData.ol.HiC.pos <- 
            switch(bindingType,
                   nearestBiDirectionalPromoters={
                       promoters(annoData.ol.HiC, upstream=0, downstream=1)
                   },
                   startSite={
                       promoters(annoData.ol.HiC, upstream=0, downstream=1)
                   }, 
                   endSite={
                       tmp <- annoData.ol.HiC
                       strand(tmp) <- ifelse(strand(tmp)=="-", "+", "-")
                       promoters(tmp, upstream=0, downstream=1)
                   },
                   stop("Not supported binding type", bindingType))
        HiC_FOR_GR.anno <- HiC_FOR_GR[subjectHits(ol)]
        HiC_REV_GR.anno <- HiC_REV_GR[subjectHits(ol)]
        annoData.ol.HiC$point_A <- start(HiC_FOR_GR.anno)
        annoData.ol.HiC$point_B <- end(HiC_FOR_GR.anno)
        annoData.ol.HiC$point_C <- start(HiC_REV_GR.anno)
        annoData.ol.HiC$point_D <- end(HiC_REV_GR.anno)
        annoData.ol.HiC$point_X <- start(annoData.ol.HiC.pos)
        annoData.ol.HiC$HiC.idx <- subjectHits(ol)
        annoData.ol.HiC
    })
    if(all(elementNROWS(anno.ol.HiCdata)==0)){
        return(GRanges())
    }
    HiC.idx.peaks <- unique(unlist(lapply(peaks.ol.HiCdata, 
                                          function(.ele) .ele$HiC.idx), 
                                   use.names = FALSE))
    HiC.idx.anno <- unique(unlist(lapply(anno.ol.HiCdata, 
                                         function(.ele) .ele$HiC.idx),
                                  use.names = FALSE))
    HiC.idx <- intersect(HiC.idx.peaks, HiC.idx.anno)
    if(length(HiC.idx)==0){
        return(GRanges())
    }
    anno.refined <- list()
    rotate.gr <- function(gr, anchor){
        strand(gr) <- ifelse(strand(gr)=="-", "+", "-")
        off.pos <- end(gr) > anchor
        start(gr[off.pos]) <- anchor[off.pos] - width(gr[off.pos])
        end(gr[off.pos]) <- anchor[off.pos]
        end(gr[!off.pos]) <- anchor[!off.pos] + width(gr[!off.pos])
        start(gr[!off.pos]) <- anchor[!off.pos]
        gr
    }
    rev.gr <- function(gr, p1, p2, px){
        tmp.shift <- p1 + p2 - 2*px
        tmp <- shift(gr, shift=tmp.shift)
        rotate.gr(tmp, px + tmp.shift)
    }
    addListInfo <- function(l, info, infoname){
        lapply(l, function(.ele){
            if(length(.ele)>0)
                mcols(.ele)[, infoname] <- rep(info, length(.ele))
            .ele
        })
    }
    unList1level <- function(l){
        l.offs <- unique(unlist(lapply(l, names)))
        sapply(l.offs, function(.ele){
            .ele <- lapply(l, `[[`, .ele)
            .ele <- .ele[sapply(.ele, length)>0]
            if(length(.ele)>0) unlist(GRangesList(.ele), use.names = FALSE)
            else NULL
        }, simplify = FALSE)
    }
    
    for(i in 1:length(anno.ol.HiCdata)){
        this.name <- names(anno.ol.HiCdata)[i]
        this.data <- anno.ol.HiCdata[[i]]
        this.data <- this.data[this.data$HiC.idx %in% HiC.idx]
        if(length(this.data)>0){
            AC.rev <- rev.gr(this.data, this.data$point_C,
                             this.data$point_A, this.data$point_X)
            AC <- shift(this.data, shift=this.data$point_C - this.data$point_A)
            CA <- shift(this.data, shift=this.data$point_A - this.data$point_C)
            AD.rev <- rev.gr(this.data, this.data$point_D,
                             this.data$point_A, this.data$point_X)
            AD <- shift(this.data, shift=this.data$point_D - this.data$point_A)
            DA <- shift(this.data, shift=this.data$point_A - this.data$point_D)
            BC.rev <- rev.gr(this.data, this.data$point_C,
                             this.data$point_B, this.data$point_X)
            BC <- shift(this.data, shift=this.data$point_C - this.data$point_B)
            CB <- shift(this.data, shift=this.data$point_B - this.data$point_C)
            BD.rev <- rev.gr(this.data, this.data$point_D,
                             this.data$point_B, this.data$point_X)
            BD <- shift(this.data, shift=this.data$point_D - this.data$point_B)
            DB <- shift(this.data, shift=this.data$point_B - this.data$point_D)
            shiftAnn <- 
                switch(this.name,
                       A.ups=list(AC=list(A.ups=NULL, 
                                          AB=AC.rev, 
                                          BC=AC.rev,
                                          CD=AC, 
                                          D.dws=AC),
                                  AD=list(A.ups=NULL, 
                                          AB=AD.rev,
                                          BC=AD.rev,
                                          CD=AD.rev,
                                          D.dws=AD),
                                  BC=list(A.ups=NULL, 
                                          AB=NULL,
                                          BC=BC.rev,
                                          CD=BC,
                                          D.dws=BC),
                                  BD=list(A.ups=NULL, 
                                          AB=NULL,
                                          BC=BD.rev,
                                          CD=BD.rev,
                                          D.dws=BD)),
                       AB=list(AC=list(A.ups=AC.rev, 
                                       AB=c(AC, CA), 
                                       BC=AC,
                                       CD=AC.rev, 
                                       D.dws=AC.rev),
                               AD=list(A.ups=AD.rev, 
                                       AB=c(AD, DA),
                                       BC=AD,
                                       CD=AD,
                                       D.dws=AD.rev),
                               BC=list(A.ups=NULL, 
                                       AB=NULL,
                                       BC=BC.rev,
                                       CD=BC,
                                       D.dws=BC),
                               BD=list(A.ups=NULL, 
                                       AB=NULL,
                                       BC=BD.rev,
                                       CD=BD.rev,
                                       D.dws=BD)),
                       BC=list(AC=list(A.ups=AC.rev, 
                                       AB=CA, 
                                       BC=c(AC, CA),
                                       CD=AC.rev, 
                                       D.dws=AC.rev),
                               AD=list(A.ups=AD.rev, 
                                       AB=DA,
                                       BC=c(AD, DA),
                                       CD=AD,
                                       D.dws=AD.rev),
                               BC=list(A.ups=BC.rev, 
                                       AB=BC.rev,
                                       BC=c(BC, CB),
                                       CD=BC.rev,
                                       D.dws=BC.rev),
                               BD=list(A.ups=BD.rev, 
                                       AB=BD.rev,
                                       BC=c(BD, DB),
                                       CD=BD,
                                       D.dws=BD.rev)),
                       CD=list(AC=list(A.ups=CA,
                                       AB=AC.rev,
                                       BC=AC.rev,
                                       CD=NULL,
                                       D.dws=NULL),
                               AD=list(A.ups=AD.rev,
                                       AB=DA,
                                       BC=DA,
                                       CD=c(AD, DA),
                                       D.dws=AD.rev),
                               BC=list(A.ups=CB,
                                       AB=CB,
                                       BC=BC.rev,
                                       CD=NULL, 
                                       D.dws=NULL),
                               BD=list(A.ups=BD.rev,
                                       AB=BD.rev,
                                       BC=DB,
                                       CD=c(BD, DB), 
                                       D.dws=BD.rev)),
                       D.dws=list(AC=list(A.ups=CA,
                                          AB=AC.rev,
                                          BC=AC.rev,
                                          CD=NULL,
                                          D.dws=NULL),
                                  AD=list(A.ups=DA,
                                          AB=AD.rev,
                                          BC=AD.rev,
                                          CD=AD.rev,
                                          D.dws=NULL),
                                  BC=list(A.ups=CB,
                                          AB=CB,
                                          BC=BC.rev,
                                          CD=NULL, 
                                          D.dws=NULL),
                                  BD=list(A.ups=DB,
                                          AB=DB,
                                          BC=BD.rev,
                                          CD=BD.rev, 
                                          D.dws=NULL)))
            shiftAnn <- 
                mapply(addListInfo, 
                       shiftAnn, names(shiftAnn), 
                       infoname="cross.link.region",
                       SIMPLIFY = FALSE)
            anno.refined[[this.name]] <- unList1level(shiftAnn)
        }
    }
    anno.refined <- 
        mapply(addListInfo, 
               anno.refined, names(anno.refined), 
               infoname="raw.annotation.region",
               SIMPLIFY = FALSE)
    anno.refined <- unList1level(anno.refined)
    anno.refined <- anno.refined[names(peaks.ol.HiCdata)]
    enhancer <- mapply(function(.anno, .peaks){
        ol.HiC.idx <- intersect(.peaks$HiC.idx, .anno$HiC.idx)
        if(length(ol.HiC.idx)==0) return(NULL)
        annot <- lapply(ol.HiC.idx, function(.HiC.id){
            .a <- .anno[.anno$HiC.idx == .HiC.id]
            .p <- .peaks[.peaks$HiC.idx == .HiC.id]
            annoPeaks(.p, .a, bindingType = bindingType,
                      bindingRegion = bindingRegion,
                      ignore.peak.strand = ignore.peak.strand)
        })
        annot <- annot[sapply(annot, length)>0]
        if(length(annot)==0) return(NULL)
        annot <- unlist(GRangesList(annot))
        colnames(mcols(annot)) <- 
            gsub("feature\\.", "feature.shift.", colnames(mcols(annot)))
        annot
    }, anno.refined, peaks.ol.HiCdata, SIMPLIFY = FALSE)
    enhancer <- enhancer[sapply(enhancer, length)>0]
    if(length(enhancer)==0) return(GRanges())
    enhancer <- mapply(function(.a, .n){
        mcols(.a)[, "peak.annotation.region"] <- .n
        .a
    }, enhancer, names(enhancer), SIMPLIFY = FALSE)
    enhancer <- unlist(GRangesList(enhancer), use.names = FALSE)
    enhancer$feature.ranges <- ranges(annoData[enhancer$feature])
    enhancer$feature.strand <- strand(annoData[enhancer$feature])
    ncols <- ncol(mcols(enhancer))
    feature.col.id <- which(colnames(mcols(enhancer))=="feature")
    mcols(enhancer) <- mcols(enhancer)[, c(1:(feature.col.id-2), feature.col.id,
                                           ncols-1, ncols, 
                                           (feature.col.id+1):(ncols-2), 
                                           feature.col.id-1)]
    enhancer$point_X <- NULL
    enhancer$HiC.idx.1 <- NULL
    colnames(mcols(enhancer)) <- 
        gsub("point_", "DNAinteractive_point_", colnames(mcols(enhancer)))
    colnames(mcols(enhancer)) <- 
        gsub("HiC", "DNAinteractive", colnames(mcols(enhancer)))
    peak.gpid <- rle(enhancer$peak.oid.to.be.deleted)
    peak.gpid$values <- seq_along(peak.gpid$values)
    peak.gpid <- inverse.rle(peak.gpid)
    enhancer <- enhancer[order(peak.gpid, enhancer$distance)]
    enhancer <- enhancer[!duplicated(paste(enhancer$feature, enhancer$peak))]
    enhancer$peak.oid.to.be.deleted <- NULL
    enhancer$DNAinteractive.ranges <- 
        ranges(DNAinteractiveData[enhancer$DNAinteractive.idx])
    enhancer$DNAinteractive.blocks <- 
        DNAinteractiveData[enhancer$DNAinteractive.idx]$blocks
    enhancer$DNAinteractive_point_A <- NULL
    enhancer$DNAinteractive_point_B <- NULL
    enhancer$DNAinteractive_point_C <- NULL
    enhancer$DNAinteractive_point_D <- NULL
    enhancer$cross.link.region <- NULL
    enhancer$raw.annotation.region <- NULL
    enhancer$peak.annotation.region <- NULL
    enhancer$DNAinteractive.idx <- NULL
    enhancer
}
