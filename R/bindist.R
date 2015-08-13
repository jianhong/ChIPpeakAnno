setClass("bindist", representation(counts="integer",
                            mids="integer",
                            halfBinSize="integer",
                            bindingType="character",
                            featureType="character"),
         validity=function(object){
             re <- TRUE
             if(is.null(object@counts)) re <- "counts is empty"
             if(is.null(object@mids)) re <- "mids is empty"
             if(is.null(object@halfBinSize)) re <- "halfBinSize is empty"
             if(length(object@halfBinSize)>1) re <- "length of halfBinSize must be 1"
             if(length(object@counts)!=length(object@mids)) 
                 re <- "the length of mids and counts are not identical"
             if(!object@bindingType %in% c("TSS", "geneEnd"))
                 re <- "the bindingType must be TSS or geneEnd"
             if(!object@featureType %in% c("transcript", "exon"))
                 re <- "the featureType must be transcript or exon"
             re
         })

setMethod("$", "bindist", function(x, name) slot(x, name))
setReplaceMethod("$", "bindist",
                 function(x, name, value){
                     slot(x, name, check = TRUE) <- value
                     x
                 })
