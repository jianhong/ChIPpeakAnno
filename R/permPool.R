#' Class \code{"permPool"}
#' 
#' An object of class \code{"permPool"} represents the possible locations to do
#' permutation test.
#' 
#' 
#' @name permPool-class
#' @rdname permPool
#' @aliases permPool permPool-class permPool-method $,permPool-method
#' $<-,permPool-method
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("permPool", grs="GRangesList", N="integer")}.
#' @slot grs object of \code{"GRangesList"} The list of binding ranges
#' @slot N vector of \code{"integer"}, permutation number for each ranges
#' @seealso \link{preparePool}, \link{peakPermTest}
#' @exportClass permPool
#' @exportMethod "$" "$<-"
#' @keywords classes

setClass("permPool", representation(grs="GRangesList",
                                    N="integer"),
         validity=function(object){
             re <- TRUE
             if(length(object@grs)!=length(object@N)) 
                 re <- "the length of grs and N are not identical"
             re
         })

setMethod("$", "permPool", function(x, name) slot(x, name))
setReplaceMethod("$", "permPool",
                 function(x, name, value){
                     slot(x, name, check = TRUE) <- value
                     x
                 })