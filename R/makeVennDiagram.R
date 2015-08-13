makeVennDiagram <- function(Peaks, NameOfPeaks, maxgap=0L, minoverlap=1L, 
                            totalTest, by=c("region", "feature", "base"), 
                            ignore.strand=TRUE, connectedPeaks=c("min", "merge", "keepAll"), 
                            method=c("hyperG", "permutation"), TxDb,
                            ...){
  ###Functions to be used
  getCountsList<-function(counts){
    CountsList <- list()
    cnt <- 1
    for(i in 1:length(counts)){
      CountsList[[i]] <- seq(cnt, length.out=counts[i])
      cnt <- cnt+counts[i]
    }
    CountsList
  }
  getVennList<-function(a, NameOfPeaks, CountsList){
    .x<-lapply(NameOfPeaks, function(.ele, cl, a){
      .y<-c()
      for(i in 1:nrow(a)){
        if(a[i,.ele]!=0) .y <- c(.y, cl[[i]])
      }
      .y
    }, CountsList, a)
    names(.x)<-NameOfPeaks
    .x
  }
  plotVenn<-function(venn_cnt, vennx, otherCounts=NULL, cat.cex = 1, cat.col = "black", 
                     cat.fontface = "plain", cat.fontfamily = "serif", ...){
    op=par(mar=c(0,0,0,0))
    on.exit(par(op))
    plot.new()
    venngrid <- venn.diagram(x=vennx, filename=NULL, cat.cex = cat.cex,
                             cat.col = cat.col, cat.fontface = cat.fontface, 
                             cat.fontfamily = cat.fontfamily, ...)
    if(grepl("^count\\.", colnames(venn_cnt)[ncol(venn_cnt)]) && connectedPeaks=="keepAll"){
        n <- which(colnames(venn_cnt)=="Counts")-1
        venn_cnt.Counts.gt.1 <- rowSums(venn_cnt[ ,1:n]) > 1
        counts <- venn_cnt[-1,"Counts"]
        counts <- counts[counts!=0]
        if(!any(duplicated(counts))){## there is no duplicated counts number
            w <- 1
            for(i in 1:length(venngrid)){
                if(inherits(venngrid[[i]], what="text") && w<nrow(venn_cnt)){
                    w <- w+1
                    cnt <- as.numeric(venngrid[[i]]$label)
                    if(!is.na(cnt) && cnt!=0){
                        j <- which(venn_cnt[,n+1]==cnt)
                        if(venn_cnt.Counts.gt.1[j]){
                            labels <- venn_cnt[j, -(1:(n+1))]
                            if(nrow(venn_cnt)==4 && sum(venn_cnt[,4]) < sum(venn_cnt[,5])){
                                ## the venn diagrame will always show small one in the right
                                labels <- rev(labels)
                            }
                            labels <- paste(as.character(labels[labels>0]), collapse="/", sep="")
                            if(labels!="") venngrid[[i]]$label <- paste(venngrid[[i]]$label, "(", labels, ")", sep="")
                        }
                    }
                }
            }
        }else{##otherwise guess the order of the number, Bugs, the map order maybe wrong
            j <- 1
            map <- switch(as.character(nrow(venn_cnt)),
                          '4'=c(3, 4, 2),
                          '8'=c(5, 7, 3, 6, 8, 4, 2),
                          '16'=c(3, 4, 2, 11, 12, 16, 8, 6, 9, 10, 14, 15, 7, 5, 13),
                          '32'=c(17, 9, 5, 3, 2, 6, 18, 19, 25, 10, 13, 21, 7, 11, 4, 8, 22, 20, 27, 26, 30, 29, 23, 15, 12, 16, 24, 28, 30, 31, 32))
            for(i in 1:length(venngrid)){
                if(inherits(venngrid[[i]], what="text") && j<length(map)){
                    while(j < length(map) && as.numeric(venngrid[[i]]$label) != venn_cnt[map[j],n+1]) j <- j+1
                    if(as.numeric(venngrid[[i]]$label) == venn_cnt[map[j],n+1] && j <= length(map)){
                        if(venn_cnt.Counts.gt.1[map[j]]){
                            labels <- venn_cnt[map[j], -(1:(n+1))]
                            labels <- paste(as.character(labels[labels>0]), collapse="/", sep="")
                            if(labels!="") venngrid[[i]]$label <- paste(venngrid[[i]]$label, "(", labels, ")", sep="")
                        }
                        j <- j+1
                    }
                }
            }
        }
    }
    if(!is.null(otherCounts)){
      tmp <- textGrob(label=otherCounts, x=0.9, y=0.1, gp=gpar(col = cat.col, cex = cat.cex, 
                                                                    fontface = cat.fontface, fontfamily = cat.fontfamily))
      venngrid <- gList(venngrid, tmp)
    }
    grid.draw(venngrid)
  }
  getPval<-function(venn_cnt, totalTest){
    n <- which(colnames(venn_cnt)=="Counts") - 1 ##ncol(venn_cnt)-1
    s <- apply(venn_cnt[,1:n], 1, sum)
    venn_cnt_s <- venn_cnt[s==2, , drop=FALSE]
	  for(i in 1:nrow(venn_cnt_s)){
		  venn_cnt_s[i,n+1] <- sum(venn_cnt[as.numeric(venn_cnt[,1:n]%*%venn_cnt_s[i,1:n])==2,n+1,drop=TRUE])
	  }
    cnt <- venn_cnt[,n+1]
    cnt_s <- apply(venn_cnt[,1:n], 2, function(.ele){
      sum(cnt[as.logical(.ele)])
    })
    p.value <- apply(venn_cnt_s, 1, function(.ele){
      ab <- cnt_s[as.logical(.ele[1:n])]
      a <- ab[1]
      b <- ab[2]
      a.and.b <- .ele[(n+1)]
      phyper <- phyper(a.and.b-1, b, totalTest-b, a, lower.tail=FALSE, log.p=FALSE)
    })
    cbind(venn_cnt_s[, 1:n, drop=FALSE], pval=p.value)
  }
  ###check inputs
  method <- match.arg(method)
  if (missing(totalTest) && method=="hyperG")
  {
    message("Missing totalTest! totalTest is required for HyperG test. 
If totalTest is missing, pvalue will be calculated by estimating 
the total binding sites of encoding region of human.
totalTest = humanGenomeSize * (2%(codingDNA) + 1%(regulationRegion)) / ( 2 * averagePeakWidth )
          = 3.3e+9 * 0.03 / ( 2 * averagePeakWidth)
          = 5e+7 /averagePeakWidth")
  }
  if (missing(Peaks))
  {
    stop("Missing Peaks which is a list of peaks in RangedData or GRanges!")
  }
  connectedPeaks <- match.arg(connectedPeaks)
  by <- match.arg(by)
  n1 = length(Peaks)
  olout_flag <- FALSE
  if(class(Peaks)=="overlappingPeaks"){
    if(!is.null(Peaks$venn_cnt) && class(Peaks$venn_cnt)=="VennCounts"){
        venn_cnt <- Peaks$venn_cnt
        n1 <- which(colnames(venn_cnt)=="Counts") - 1 ##ncol(venn_cnt)-1
        if(missing(NameOfPeaks)){
            NameOfPeaks <- colnames(venn_cnt)[1:n1]
        }
        olout_flag <- TRUE
    }else{
        stop("Input is object of overlappingPeaks. But it does not have venn counts data.")
    }
  }
  if (missing(NameOfPeaks) ||  mode(NameOfPeaks) != "character")
  {
    warning("Missing required character vector NameOfPeaks. NameOfPeaks will be extract from the names of input Peaks.")
    NameOfPeaks <- names(Peaks)
    if(is.null(NameOfPeaks)){
      NameOfPeaks <- paste("peaks", 1:length(Peaks), sep="")
    }
  }
  n2 = length(NameOfPeaks)
  if (n1 <2)
  {
    stop("The number of element in Peaks is less than 1, need at least 2 peak ranges!")
  }
  if (n1 > 5) stop("The length of input peaks list should no more than 5")
  if (n1 > n2)
  {
     stop("The number of element in NameOfPeaks is less than the number of elements in Peaks, need to be equal!")
  }
  if (n1 < n2)
  {
    warning("The number of element in NameOfPeaks is larger than the number of elements in Peaks! NameOfPeaks will be cut by the length of Peaks")
    NameOfPeaks <- NameOfPeaks[1:n1]
  }
  NameOfPeaks <- make.names(NameOfPeaks, unique=TRUE, allow_=TRUE)
  if(!olout_flag){
      for(i in seq_len(n1)){
        if (!inherits(Peaks[[i]], c("RangedData", "GRanges")))
        {
          err.msg = paste("Element", i, "in Peaks is not a valid RangedData or GRanges object", sep=" ")
          stop(err.msg)
        }
        if (inherits(Peaks[[i]], "RangedData")) Peaks[[i]] <- toGRanges(Peaks[[i]], format="RangedData")
        if (is.null(names(Peaks[[i]]))) names(Peaks[[i]]) <-  formatC(1:length(Peaks[[i]]), width=nchar(length(Peaks[[i]])), flag="0")
      }
      if(!missing(totalTest)){
          if(totalTest < max(sapply(Peaks, length))){
              stop("totaltest specifies the total number of possible peaks in the testing space. It should be larger than the largest peak numbers in the input sets. Please see more details at http://pgfe.umassmed.edu/ChIPpeakAnno/FAQ.html")
          }
      }
      names(Peaks) <- NameOfPeaks
      venn_cnt <- getVennCounts(Peaks, maxgap=maxgap, minoverlap=minoverlap, by=by, ignore.strand=ignore.strand, connectedPeaks=connectedPeaks)
  }
  colnames(venn_cnt)[1:n1] <- NameOfPeaks 
  Counts <- getCountsList(venn_cnt[,"Counts"])
  vennx <- getVennList(venn_cnt, NameOfPeaks, Counts)
  if(method=="hyperG"){
      if(!missing(totalTest)){
          otherCount <- totalTest - sum(venn_cnt[, "Counts"])
          venn_cnt[1, "Counts"] <- otherCount
          p.value <- getPval(venn_cnt, totalTest)
      }else{
          otherCount <- venn_cnt[1, "Counts"]
          if(class(Peaks)=="overlappingPeaks"){
              averagePeakWidth = median(unlist(lapply(Peaks$peaklist, width)))
          }else{
              averagePeakWidth = median(unlist(lapply(Peaks, width))) 
          }
          
          p.value <- getPval(venn_cnt, round(5e+7/averagePeakWidth))
      }
  }else{
      ##method=="permutation"
      p.value <- NA
      if(olout_flag){
          warning("Input is an object of overlappingPeaks. Can not do permutation test. Please try ?peakPermTest")
      }else{
          if(missing(TxDb)){
              warning("TxDb is missing. Please try ?peakPermTest later.")
          }else{
              n <- which(colnames(venn_cnt)=="Counts") - 1 ##ncol(venn_cnt)-1
              s <- apply(venn_cnt[,1:n], 1, sum)
              venn_cnt_s <- venn_cnt[s==2, , drop=FALSE]
              p.value <- apply(venn_cnt_s, 1, function(.ele){
                  AB <- Peaks[as.logical(.ele)]
                  pt <- peakPermTest(AB[[1]], AB[[2]], TxDb=TxDb, ...)
                  pt$pval
              })
              p.value <- cbind(venn_cnt_s[, 1:n, drop=FALSE], pval=p.value)
          }
      }
      if(!missing(totalTest)){
          otherCount <- totalTest - sum(venn_cnt[, "Counts"])
          venn_cnt[1, "Counts"] <- otherCount
      }else{
          otherCount <- venn_cnt[1, "Counts"]
      }
  }
  plotVenn(venn_cnt, vennx, otherCount=otherCount, ...)
  return(list(p.value=p.value, vennCounts=venn_cnt))
}