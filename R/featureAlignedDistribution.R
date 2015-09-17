featureAlignedDistribution <- function(cvglists, feature.gr, 
                                       upstream, downstream, 
                                       n.tile=100, zeroAt, ...){
    stopifnot(class(feature.gr)=="GRanges")
    dots <- list(...)
    
    grWidr <- unique(width(feature.gr))
    if(missing(upstream) || missing(downstream)){
        if(length(grWidr)!=1){
            stop("The width of feature.gr is not identical.")
        }
        if(missing(zeroAt)) {
            zeroAt <- 0.5
            message("zero is set as the center of the feature.gr")
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
        zeroAt <- upstream/(upstream + downstream)
        end(feature.gr) <- start(feature.gr) + downstream
        start(feature.gr) <- start(feature.gr) - upstream
        grWidr <- unique(width(feature.gr))
    }
    stopifnot(is.numeric(zeroAt))
    stopifnot(zeroAt>=0)
    if(zeroAt<=1){
        zero <- round(grWidr*zeroAt)
    }else{
        zero <- round(zeroAt)
    }
    
    grWid <- c(0, grWidr) - zero
    grWidLab <- grid.pretty(grWid)
    grWidAt <- (grWidLab+zero)/grWidr*n.tile
    if(inherits(cvglists, c("SimpleRleList", "RleList", "CompressedRleList"))){
        cvglistsName <- substitute(deparse(cvglists))
        cvglists <- list(cvglists)
        names(cvglists) <- cvglistsName
    }
    if(class(cvglists)!="list"){
        stop("cvglists must be output of featureAlignedSignal or 
             a list of SimpleRleList or RleList")
    }
    cls <- sapply(cvglists, class)
    if(all(cls=="matrix")){
        cov <- cvglists
        if(ncol(cov[[1]])!=n.tile){
            stop("n.tile must keep same as featureAligendSignal.")
        }
    }else{
        if(any(!cls %in% c("SimpleRleList", "RleList", "CompressedRleList")))
            stop("cvglists must be a list of SimpleRleList or RleList")
        cov <- featureAlignedSignal(cvglists, feature.gr, n.tile=n.tile)
    }
    ## normalized read density
    density <- sapply(cov, colMeans)
    matplot(density, ..., xaxt="n")
    axis(1, at = grWidAt, labels = grWidLab)
    lty <- if(!is.null(dots$lty)) dots$lty else 1:5
    lwd <- if(!is.null(dots$lwd)) dots$lwd else 1
    col <- if(!is.null(dots$col)) dots$col else 1:6
    legend("topright", legend=colnames(density), col=col,
           lty=lty, lwd=lwd)
    return(invisible(density))
}