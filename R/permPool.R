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