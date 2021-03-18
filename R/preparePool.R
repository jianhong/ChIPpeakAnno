#' prepare data for permutation test
#' 
#' prepare data for permutation test \link{peakPermTest}
#' 
#' 
#' @param TxDb an object of \link[GenomicFeatures:TxDb-class]{TxDb}
#' @param template an object of \link[GenomicRanges:GRanges-class]{GRanges}
#' @param bindingDistribution an object of \link{bindist}
#' @param bindingType the relevant position to features
#' @param featureType feature type, transcript or exon.
#' @param seqn seqnames. If given, the pool for permutation will be restrict in
#' the given chromosomes.
#' @return a list with two elements, grs, a list of
#' \link[GenomicRanges:GRanges-class]{GRanges}.  N, the numbers of elements
#' should be drawn from in each GRanges.
#' @author Jianhong Ou
#' @seealso \link{peakPermTest}, \link{bindist}
#' @keywords misc
#' @export
#' @importFrom S4Vectors queryHits subjectHits
#' @examples
#' 
#'     if(interactive() || Sys.getenv("USER")=="jianhongou"){
#'         path <- system.file("extdata", package="ChIPpeakAnno")
#'         peaksA <- toGRanges(file.path(path, "peaks.narrowPeak"), 
#'                             format="narrowPeak")
#'         peaksB <- toGRanges(file.path(path, "MACS2_peaks.xls"), format="MACS2")
#'         library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#'         ppp <- preparePool(TxDb.Hsapiens.UCSC.hg19.knownGene, 
#'                            peaksA, bindingType="TSS",
#'                            featureType="transcript")
#'     }
#' 
preparePool <- function(TxDb, template, bindingDistribution,
                        bindingType=c("TSS", "geneEnd"), 
                        featureType=c("transcript", "exon"), seqn=NA){
    if(missing(TxDb)){
        stop("TxDb is requried!")
    }
    if(!is(TxDb, "TxDb")){
        stop("TxDb must be an object of TxDb")
    }
    if(!missing(bindingDistribution)){
        if(!is(bindingDistribution, "bindist")){
            stop("bindingDistribution must be an object of bindist")
        }
    }
    bindingType <- match.arg(bindingType)
    featureType <- match.arg(featureType)
    if(featureType=="transcript"){
        tx <- transcripts(x=TxDb)
        tx <- unique(tx)
        names(tx) <- seq_along(tx)
    }else{
        tx <- exons(x=TxDb)
        tx <- unique(tx)
        names(tx) <- seq_along(tx)
    }
    
    tx <- unique(tx)
    
    if(missing(bindingDistribution)){
        bindingDistribution <- 
            buildBindingDistribution(template, AnnotationData=tx,
                                     bindingType=bindingType,
                                     featureType=featureType)
    }
    offset <- bindingDistribution$mids
    diff <- bindingDistribution$halfBinSize
    N <- bindingDistribution$counts
    
    if(!is.na(seqn[1])) tx <- tx[as.character(seqnames(tx)) %in% seqn]
    
    grs <- lapply(offset, function(.off){
        x <- NULL
        suppressWarnings({
            x <- shift(tx, .off)
            start(x) <- start(x) - diff
            width(x) <- 2*diff
            x <- trim(x)
            x <- x[width(x) > 1]
        })
        ol <- findOverlaps(x, tx)
        if(length(ol)>0){
            if(.off+diff<0){
                x <- x[-(unique(queryHits(ol)))]
            }else{
                x <- x[-(unique(queryHits(ol)[
                    names(x)[queryHits(ol)]!=names(tx)[subjectHits(ol)]]))]
            }
        }
    })
    
    list(grs=GRangesList(grs), N=as.integer(N))
}
