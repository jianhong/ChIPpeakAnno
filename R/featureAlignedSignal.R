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
    feature.gr <- promoters(feature.gr, upstream = upstream,
                            downstream = downstream + 1)
    grWidr <- unique(width(feature.gr))
  }
  if(any(start(feature.gr)<1)){
    warning("Some start position of the peaks are less than 1!",
            "They will be filtered.")
    feature.gr <- feature.gr[start(feature.gr)>0]
  }
  if(length(feature.gr)<2){
    stop("Length of feature.gr less than 2.")
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
  seqLen <-  lapply(cvglists, function(.ele) 
    sapply(.ele, function(.e) sum(runLength(.e))))
  seqLen.keep <- table(unlist(sapply(seqLen, names)))==length(cvglists)
  seqLen <- seqLen[[1]][seqLen.keep]
  seqLen <- seqLen[!is.na(seqLen)]
  if(length(seqLen)>0){
    feature.gr.subset <- 
      feature.gr[seqnames(feature.gr) %in% names(seqLen)]
    if(any(end(feature.gr.subset)>seqLen[as.character(seqnames(feature.gr.subset))])){
      warning("Some end position of the peaks are out of bound!",
              "They will be filtered.")
      feature.gr <- 
        feature.gr.subset[end(feature.gr.subset) <=
                            seqLen[as.character(seqnames(feature.gr.subset))]]
    }
  }
  #feature.gr.bck <- feature.gr
  stopifnot(is.numeric(n.tile))
  n.tile <- round(n.tile)
  grL <- tile(feature.gr, n=n.tile)
  idx <- as.character(strand(feature.gr))=="-"
  if(length(idx)>0) {
      grL[idx] <- as(lapply(grL[idx], rev), "GRangesList")
  }
  grL.len <- lengths(grL)
  grL <- unlist(grL)
  grL$oid <- rep(seq_along(feature.gr), grL.len)
  grL$nid <- unlist(lapply(grL.len, seq.int))
  grL.len <- length(grL)
  grL.s <- split(grL, as.character(seqnames(grL)))
  grL <- unlist(grL.s) ## make sure grL and grL.s keep the same order
  seqn <- names(grL.s)
  seql <- seqlengths(feature.gr)
  seql <- seql[seqn]
  seql.f <- range(feature.gr)
  seql.f <- seql.f[match(seqn, seqnames(seql.f))]
  seql[is.na(seql)] <- width(seql.f)[is.na(seql)] + 1
  trimChar <- function(x, width){
      len <- nchar(x)
      if(len <= width) return(x)
      return(paste0(strtrim(x, width=width), "..."))
  }
  rowname.feature.gr <- paste0(as.character(seqnames(feature.gr)), ":",
                               start(feature.gr), "-",
                               end(feature.gr))
  cov <- lapply(cvglists, function(.dat){
      .dat <- .dat[seqn[seqn %in% names(.dat)]]
      if(length(.dat)!=length(seqn)){
          warning(paste(seqn[!seqn %in% names(.dat)], collapse=", "), 
                  ifelse(length(seqn[!seqn %in% names(.dat)])>1, " are", " is"), 
                  " not in cvglists. seqlevels of cvglist are ", 
                  trimChar(paste(names(.dat), collapse=", "), width=60))
          for(i in seqn[!seqn %in% names(.dat)]){
              .dat[[i]] <- Rle(0, seql[i])
          }
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
      .dat <- .dat[seqn]
      vw <- Views(as(.dat, "RleList"), grL.s)
      vm <- viewMeans(vw)
      vm <- unlist(vm)
      stopifnot(length(vm)==grL.len)
      mm <- matrix(0, nrow=length(feature.gr), ncol=n.tile)
      mm[grL$oid+length(feature.gr)*(grL$nid-1)] <- vm
      rownames(mm) <- rowname.feature.gr
      mm
  })
  
  cov
}