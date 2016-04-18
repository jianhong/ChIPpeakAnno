featureAlignedSignal <- function(cvglists, feature.gr, 
                                 upstream, downstream, 
                                 n.tile=100, ...){
    stopifnot(class(feature.gr)=="GRanges")
    
    grWidr <- unique(width(feature.gr))
    if(missing(upstream) || missing(downstream)){
        if(length(grWidr)!=1){
            stop("The width of feature.gr is not identical.")
        }
    }else{
        if(!is.numeric(upstream) || !is.numeric(downstream)){
            stop("upstream and downstream must be integers")
        }
        if(upstream<0 || downstream<0){
            stop("upstream and downstream must be not less than 0")
        }
        upstream <- as.integer(upstream)
        downstream <- as.integer(downstream)
        if(length(grWidr)!=1){
            start(feature.gr) <- start(feature.gr)+floor(grWidr/2)
            width(feature.gr) <- 1
            warning("feature.gr is set to the center of feature.gr")
        }
        end(feature.gr) <- start(feature.gr) + downstream
        start(feature.gr) <- start(feature.gr) - upstream
        grWidr <- unique(width(feature.gr))
    }
    if(inherits(cvglists, c("SimpleRleList", "RleList", "CompressedRleList"))){
        cvglistsName <- substitute(deparse(cvglists))
        cvglists <- list(cvglists)
        names(cvglists) <- cvglistsName
    }
    if(class(cvglists)!="list"){
        stop("cvglists must be a list of SimpleRleList or RleList")
    }
    cls <- sapply(cvglists, class)
    if(any(!cls %in% c("SimpleRleList", "RleList", "CompressedRleList")))
        stop("cvglists must be a list of SimpleRleList or RleList")
    feature.gr.bck <- feature.gr
    feature.gr$oid <- 1:length(feature.gr)
    feature.gr <- split(feature.gr, as.character(seqnames(feature.gr)))
    stopifnot(is.numeric(n.tile))
    n.tile <- round(n.tile)
    grL <- lapply(feature.gr, tile, n=n.tile)
    oid <- lapply(feature.gr, function(.ele) .ele$oid)
    oid <- unlist(oid)
    
    cov <- lapply(cvglists, function(.dat){
        do.call(rbind, 
                mapply(function(x, n){
                    .d <- .dat[[n]]
                    if(length(.d)==0){
                        message(n, " is not in cvglists. seqlevels of cvglist is ", names(.dat))
                        return(NULL)
                    }
                    do.call(rbind, lapply(x, function(.ele){
                        log2(viewMeans(Views(.d, ranges(.ele)))+1)
                    }))
                }, grL, names(grL), SIMPLIFY=FALSE))
    })
    lapply(cov, function(.ele) .ele[order(oid), ])
}