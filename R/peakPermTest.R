peakPermTest <- function(peaks1, peaks2, ntimes=100, 
                         seed=as.integer(Sys.time()),
                         mc.cores=getOption("mc.cores", 2L),
                         maxgap=-1L, pool,
                         TxDb, bindingDistribution,
                         bindingType=c("TSS", "geneEnd"), 
                         featureType=c("transcript", "exon"),
                         seqn=NA, ...){
    if(!is(peaks1, "GRanges") || !is(peaks2, "GRanges")){
        stop("class of peaks1 and peaks2 must be GRanges.")
    }
    if(!missing(bindingDistribution)){
        if(class(bindingDistribution)!="bindist"){
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
        if(class(pool)!="permPool"){
            stop("class of pool must be permPool")
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
