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
        if(length(grWidr)!=1 || any(grWidr!=1)){
            start(feature.gr) <- start(feature.gr)+floor(width(feature.gr)/2)
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
    #feature.gr.bck <- feature.gr
    feature.gr$oid <- 1:length(feature.gr)
    feature.gr <- split(feature.gr, as.character(seqnames(feature.gr)))
    stopifnot(is.numeric(n.tile))
    n.tile <- round(n.tile)
    grL <- lapply(feature.gr, function(.ele){
        .gr <- tile(.ele, n=n.tile)
        .gr <- unlist(.gr)
        .gr$oid <- rep(.ele$oid, each=n.tile)
        .gr
        })
    seqn <- names(grL)
    cov <- lapply(cvglists, function(.dat){
        .dat <- .dat[seqn]
        .dat.eleLen <- elementNROWS(.dat)
        if(any(.dat.eleLen==0)){
            warning(paste(names(.dat)[.dat.eleLen], collapse=", "), 
                    " is not in cvglists. seqlevels of cvglist is ", 
                    names(.dat))
            .dat <- .dat[.dat.eleLen!=0]
        }
        warn <- sapply(.dat, function(.ele){
            any(is.na(runValue(.ele)))
        })
        if(any(warn)){
            warning("cvglists contain NA values.")
        }
        .dat <- sapply(.dat, function(.ele){
            if(any(is.infinite(runValue(.ele)))){
                warning("cvglists contain infinite values. ", 
                        "infinite value will be converted to NA.")
                runValue(.ele)[is.infinite(runValue(.ele))] <- NA
            }
            .ele
        })
        .viewmean <- mapply(function(.d, .gr){
            .m <- matrix(viewMeans(Views(.d, ranges(.gr)), na.rm=TRUE), 
                   ncol=n.tile, byrow = TRUE)
            rownames(.m) <- as.character(matrix(.gr$oid, ncol=n.tile, 
                                                byrow=TRUE)[, 1])
            .m
        }, .dat, grL[names(.dat)], SIMPLIFY=FALSE)
        do.call(rbind, .viewmean)
    })
    lapply(cov, function(.ele) unname(.ele[order(as.numeric(rownames(.ele))), ]))
}