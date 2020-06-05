#' plot distribution in given ranges
#' 
#' plot distribution in the given feature ranges
#' 
#' 
#' @param cvglists Output of \link{featureAlignedSignal} or a list of
#' \link[IRanges:AtomicList-class]{SimpleRleList} or
#' \link[IRanges:AtomicList-class]{RleList}
#' @param feature.gr An object of \link[GenomicRanges:GRanges-class]{GRanges}
#' with identical width.  If the width equal to 1, you can use upstream and
#' downstream to set the range for plot.  If the width not equal to 1, you can
#' use zeroAt to set the zero point of the heatmap.
#' @param upstream,downstream upstream or dwonstream from the feature.gr.
#' @param zeroAt zero point position of feature.gr
#' @param n.tile The number of tiles to generate for each element of
#' feature.gr, default is 100
#' @param ... any paramters could be used by \link[graphics]{matplot}
#' @return invisible matrix of the plot.
#' @author Jianhong Ou
#' @seealso See Also as \link{featureAlignedSignal},
#' \link{featureAlignedHeatmap}
#' @keywords misc
#' @export
#' @importFrom BiocGenerics start end width strand
#' @importFrom graphics matplot
#' @examples
#' 
#'   cvglists <- list(A=RleList(chr1=Rle(sample.int(5000, 100), 
#'                                       sample.int(300, 100))), 
#'                    B=RleList(chr1=Rle(sample.int(5000, 100), 
#'                                       sample.int(300, 100))))
#'   feature.gr <- GRanges("chr1", IRanges(seq(1, 4900, 100), width=100))
#'   featureAlignedDistribution(cvglists, feature.gr, zeroAt=50, type="l")
#' 
featureAlignedDistribution <- function(cvglists, feature.gr, 
                                       upstream, downstream, 
                                       n.tile=100, zeroAt, ...){
    stopifnot(is(feature.gr, "GRanges"))
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
        if(length(grWidr)!=1 || any(grWidr!=1)){
            start(feature.gr) <- start(feature.gr)+floor(width(feature.gr)/2)
            width(feature.gr) <- 1
            warning("feature.gr is set to the center of feature.gr")
        }
        if(!missing(zeroAt)){
            warning("zeroAt will be ignored.")
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
    if(!is(cvglists, "list")){
        stop("cvglists must be output of featureAlignedSignal or 
             a list of SimpleRleList or RleList")
    }
    
    cls <- sapply(cvglists, is, class2="matrix")
    if(all(cls)){
      cov <- cvglists
      if(ncol(cov[[1]])!=n.tile){
        stop("n.tile must keep same as featureAligendSignal.")
      }
    }else{
      cls <- sapply(cvglists, inherits, 
                    what = c("SimpleRleList", "RleList", "CompressedRleList"))
      if(any(!cls))
        stop("cvglists must be a list of SimpleRleList or RleList")
      cov <- featureAlignedSignal(cvglists, feature.gr, n.tile=n.tile)
    }
    
    ## normalized read density
    if(any(sapply(cov, function(.ele) any(is.na(.ele))))){
        warning("cvglists contain NA values. ", 
                "NA value will be omit.")
    }
    density <- sapply(cov, colMeans, na.rm=TRUE)
    try({
        matplot(density, ..., xaxt="n")
        axis(1, at = grWidAt, labels = grWidLab)
        lty <- if(!is.null(dots$lty)) dots$lty else 1:5
        lwd <- if(!is.null(dots$lwd)) dots$lwd else 1
        col <- if(!is.null(dots$col)) dots$col else 1:6
        legend("topright", legend=colnames(density), col=col,
               lty=lty, lwd=lwd)
    })
    
    return(invisible(density))
}
