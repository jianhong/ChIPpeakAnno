getAllPeakSequence <- function(myPeakList, 
                               upstream=200L, downstream=upstream, 
                               genome, AnnotationData)
{
    if (!inherits(myPeakList, c("RangedData", "GRanges"))) {
        stop("No valid myPeakList passed in. It needs to be GRanges object")
    }
    if (missing(genome))
    {
        stop("genome is required parameter, 
             please pass in either a BSgenome object or a Mart object!")
    }
    myPeakList.bk <- myPeakList
    if(inherits(myPeakList, "RangedData")){
        myPeakList <- RangedData2GRanges(myPeakList)
    }
    if (is(genome, "BSgenome"))
    {
        strand <- strand(myPeakList)
        width <- width(myPeakList)
        ##  start(myPeakList)
        start(myPeakList) <- ifelse(strand=="-", 
                                    start(myPeakList) - as.numeric(downstream),
                                    start(myPeakList) - as.numeric(upstream))
        lt0 <- start(myPeakList)<=0
        start <- start(myPeakList)-1
        start[!lt0] <- 0
        start(myPeakList)[lt0] <- 1
        width(myPeakList) <- width + as.numeric(upstream) + 
            as.numeric(downstream) + start
        strand[strand!="-"] <- "+"
        strand(myPeakList) <- strand
        ##how to make the seqname be same is a question.
        if(!all(seqlevels(myPeakList) %in% seqnames(genome))){
            myPeakList <- formatSeqnames(myPeakList)
            genome <- formatSeqnames(genome)
        }
        ends <- unlist(
            apply(
                cbind(end(myPeakList), 
                      seqlengths(genome)[as.character(seqnames(myPeakList))]), 
                1,
                min))
        ends <- unname(ifelse(is.na(ends), end(myPeakList), ends))
        keep <- start(myPeakList) <= ends
        end(myPeakList)[keep] <- ends[keep]
        myPeakList <- myPeakList[keep]
        seq <- getSeq(genome, myPeakList, as.character=TRUE)
        
        myPeakList <- myPeakList.bk
        myPeakList$upstream <- rep(upstream,length(myPeakList.bk))
        myPeakList$downstream <- rep(downstream,length(myPeakList.bk))
        myPeakList$sequence <- seq
        if(all(seqlevels(myPeakList) %in% seqlevels(genome))){
            seqlengths(myPeakList) <- seqlengths(genome)[seqlevels(myPeakList)]
        }   

        myPeakList
    }else if (is(genome, "Mart")){
        if (missing(AnnotationData)) {
            message("No AnnotationData as GRanges is passed in, 
                    so now querying biomart database for AnnotationData ....")
            AnnotationData <- getAnnotation(genome)
            message("Done querying biomart database, start annotating ....
                    Better way would be calling getAnnotation before 
                    querying for sequence")
        }
        if (!inherits(AnnotationData, c("RangedData", "GRanges"))) {
            stop("AnnotationData needs to be GRanges object. 
                 Better way would be calling getAnnotation.")
        }
        if (inherits(AnnotationData, "RangedData")) {
            AnnotationData <- RangedData2GRanges(AnnotationData)
        }
        
        downstream.bk = downstream
        plusAnno = AnnotationData[strand(AnnotationData)=="+"]
        temp = annotatePeakInBatch(myPeakList, AnnotationData=plusAnno)
        TSSlength =temp$end_position - temp$start_position
        downstream = end(temp) - start(temp) + downstream
        temp$downstream=  downstream
        temp$TSSlength = TSSlength
        
        myList3 = as.data.frame(temp)
        rm(temp)
        
        if (dim(myList3)[1] != 0)
        {
            l3 =  cbind(as.character(myList3$feature), 
                        as.numeric(as.character(myList3$distancetoFeature)), 
                        as.numeric(rep(upstream, dim(myList3)[1])), 
                        as.numeric(as.character(myList3$downstream)), 
                        as.numeric(as.character(myList3$start_position)), 
                        as.numeric(as.character(myList3$end_position)))
            r3 = apply(l3, 1, getGeneSeq, genome)
        }
        else
        {
            r3 = 0
        }
        
        if (is.list(r3))
        {
            r = as.data.frame(do.call("rbind",r3))
        }
        else
        {
            stop("No sequence found error!")
        }
        colnames(r)= c("feature", "distancetoFeature", "upstream", 
                       "downstream", "seq")
        r4 = merge(r, myList3)
        GRanges(seqnames=as.character(r4$space),
                ranges=IRanges(start=r4$start, 
                               end  =r4$end, 
                               names=as.character(r4$name)),
                strand="+",
                upstream=rep(upstream, dim(r4)[1]), 
                downstream= rep(downstream.bk, dim(r4)[1]),
                sequence=unlist(r4$seq))
    }else{
        stop("genome needs to be either a BSgenome object or Mart object!")
    }
}
