#' Return the value from a Bimap objects
#' 
#' Search by name for an Bimap object.
#' 
#' 
#' @param x,envir,mode,ifnotfound,inherits see
#' \link[AnnotationDbi:Bimap-envirAPI]{mget}
#' @param output return the all or first item for each query
#' @return a character vector
#' @author Jianhong Ou
#' @seealso See Also as \code{\link[base:get]{mget}},
#' \code{\link[AnnotationDbi:Bimap-envirAPI]{mget}}
#' @keywords misc
#' @export
#' @examples
#' 
#'     library(org.Hs.eg.db)
#'     xget(as.character(1:10), org.Hs.egSYMBOL)
#' 
xget <- function(x, envir, mode, ifnotfound=NA, inherits, 
                 output=c("all", "first", "last")){
    output <- match.arg(output)
    y <- mget(x=x, envir=envir, mode=mode, 
              ifnotfound=ifnotfound, inherits=inherits)
    switch(output, 
           all=sapply(y, base::paste, collapse=";"),
           first=sapply(y, `[`, 1),
           last=sapply(y, function(.ele) .ele[length(.ele)])
    )
}
