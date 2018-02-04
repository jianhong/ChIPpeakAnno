expandGR <- function(gr, size){
    if(size == 0){
        return(gr)
    }
    if(!is(gr, "GRanges")){
        stop("gr must be an object of GRanges")
    }
    if(size < 0){
        gr <- gr[width(gr)>2*abs(size)]
    }
    
    suppressWarnings({
        start(gr) <- start(gr) - size
        end(gr) <- end(gr) + size
    })
    
    gr <- trim(gr)
    gr
}