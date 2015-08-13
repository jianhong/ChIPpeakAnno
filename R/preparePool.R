preparePool <- function(TxDb, template, bindingDistribution,
                        bindingType=c("TSS", "geneEnd"), 
                        featureType=c("transcript", "exon"), seqn=NA){
    if(missing(TxDb)){
        stop("TxDb is requried!")
    }
    if(class(TxDb)!="TxDb"){
        stop("TxDb must be an object of TxDb")
    }
    if(!missing(bindingDistribution)){
        if(class(bindingDistribution)!="bindist"){
            stop("bindingDistribution must be an object of bindist")
        }
    }
    bindingType <- match.arg(bindingType)
    featureType <- match.arg(featureType)
    if(featureType=="transcript"){
        tx <- transcripts(x=TxDb)
        tx <- unique(tx)
        names(tx) <- 1:length(tx)
    }else{
        tx <- exons(x=TxDb)
        tx <- unique(tx)
        names(tx) <- 1:length(tx)
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
    
    if(!is.na(seqn)) tx <- tx[as.character(seqnames(tx)) %in% seqn]
    
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