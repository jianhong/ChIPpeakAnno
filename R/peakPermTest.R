#' Permutation Test for two given peak lists
#' 
#' Performs a permutation test to seee if there is an association between two
#' given peak lists.
#' 
#' 
#' @param peaks1,peaks2 an object of
#' \link[GenomicRanges:GRanges-class]{GRanges}
#' @param ntimes number of permutations
#' @param seed random seed
#' @param mc.cores The number of cores to use. see mclapply.
#' @param maxgap See \link[IRanges:findOverlaps-methods]{findOverlaps} in the
#' IRanges package for a description of these arguments.
#' @param pool an object of \link{permPool}
#' @param TxDb an object of \link[GenomicFeatures:TxDb-class]{TxDb}
#' @param bindingDistribution an object of \link{bindist}
#' @param bindingType where the peaks should bind, TSS or geneEnd
#' @param featureType what annotation type should be used for detecting the
#' binding distribution.
#' @param seqn default is NA, which means not filter the universe pool for
#' sampling.  Otherwise the universe pool will be filtered by the seqnames in
#' seqn.
#' @param \dots further arguments to be passed to \link[regioneR]{numOverlaps}.
#' @return A list of class permTestResults. See \link[regioneR]{permTest}
#' @author Jianhong Ou
#' @seealso \link{preparePool}, \link{bindist}
#' @references Davison, A. C. and Hinkley, D. V. (1997) Bootstrap methods and
#' their application, Cambridge University Press, United Kingdom, 156-160
#' @keywords misc
#' @export
#' @importFrom regioneR permTest numOverlaps
#' @examples
#' 
#'     path <- system.file("extdata", package="ChIPpeakAnno")
#'     #files <- dir(path, pattern="[12]_WS170.bed", full.names=TRUE)
#'     #peaks1 <- toGRanges(files[1], skip=5)
#'     #peaks2 <- toGRanges(files[2], skip=5)
#'     #peakPermTest(peaks1, peaks2, TxDb=TxDb.Celegans.UCSC.ce6.ensGene)
#'     if(interactive()){
#'         peaks1 <- toGRanges(file.path(path, "MACS2_peaks.xls"), 
#'                             format="MACS2")
#'         peaks2 <- toGRanges(file.path(path, "peaks.narrowPeak"), 
#'                             format="narrowPeak")
#'         library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#'         peakPermTest(peaks1, peaks2, 
#'                TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene, min.pctA=10)
#'     }
#' 
peakPermTest <- function(peaks1, peaks2, ntimes=100, 
                         seed=as.integer(Sys.time()),
                         mc.cores=getOption("mc.cores", 2L),
                         maxgap=-1L, pool,
                         TxDb, bindingDistribution,
                         bindingType=c("TSS", "geneEnd"), 
                         featureType=c("transcript", "exon"),
                         seqn=NA, ...){
    if(!is(peaks1, "GRanges") || !is(peaks2, "GRanges")){
        stop("peaks1 and peaks2 must be GRanges objects.")
    }
    if(!missing(bindingDistribution)){
        if(!is(bindingDistribution, "bindist")){
            stop("bindingDistribution must be an object of bindist")
        }
    }
    bindingType <- match.arg(bindingType)
    featureType <- match.arg(featureType)
    if(missing(pool)){
        pool <- preparePool(TxDb, template=peaks1, 
                            bindingDistribution=bindingDistribution, 
                            bindingType=bindingType, 
                            featureType=featureType, seqn=seqn)
    }else{
        if(!is(pool, "permPool")){
            stop("pool must be a permPool object")
        }
    }
    set.seed(seed)
    pt <- permTest(A=peaks1, ntimes=ntimes,
                   randomize.function=randPeaks, 
                   grs=pool$grs, N=pool$N, 
                   maxgap=maxgap,
                   evaluate.function=cntOverlaps, 
                   B=peaks2, verbose=FALSE, alternative="greater", 
                   mc.set.seed=FALSE, mc.cores=mc.cores,
                   ...)
    return(pt)
}
