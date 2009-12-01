getAllPeakSequence <-
function(myPeakList, upstream=200, downstream=200, genome, AnnotationData)
{
		if (missing(genome))
		{
			stop("genome is required parameter, please pass in either a BSgenome object or a Mart object!")
		}
		if (class(genome)== "BSgenome")
		{
			Start = start(myPeakList) - as.numeric(upstream)
			End = end(myPeakList) + as.numeric(downstream)
		
			chr = as.character(space(myPeakList))
			if (!length(i<-grep("chr",chr[1])) &  length(i<- grep("chr",names(seqlengths(genome))[1])))
			{
				chr = paste("chr", chr, sep="")
			}
					
			for (i in 1:length(Start))		
			{
				thisChr =chr[i]
				thisEnd = min(End[i], seqlengths(genome)[thisChr][[1]])
				thisStart = max(1, Start[i])
				if (i ==1)
				{
					seq = getSeq(genome, thisChr, start=thisStart, end=thisEnd, width=NA, as.character=TRUE)
				}
				else
				{
					seq = c(seq, getSeq(genome, thisChr, start=thisStart, end=thisEnd, width=NA, as.character=TRUE))
				}
			}
			RangedData(IRanges(start=start(myPeakList), end = end(myPeakList), names=rownames(myPeakList)), space=chr, 
				upstream=rep(upstream,length(start(myPeakList))), downstream=rep(downstream,length(start(myPeakList))), sequence=seq)
		}
		else if (class(genome) =="Mart")
		{
			if (missing(AnnotationData))
			{
				message("No AnnotationData is passed in, so now querying biomart database for AnnotationData!")
				AnnotationData <- getAnnotation(genome)
				message("Done querying biomart database for AnnotationData, better way would be to call getAnnotation first, start querying for sequence ....")
			}
			if (class(AnnotationData)  != "RangedData")
			{
				message("No AnnotationData as RangedData is passed in, so now querying biomart database for AnnotationData!")
				AnnotationData = getAnnotation(genome)
				message("Done querying biomart database for AnnotationData, better way would be to call getAnnotation first, start querying for sequence ....")
			}
			downstream.bk = downstream
			plusAnno = AnnotationData[AnnotationData$strand==1,]
			temp = annotatePeakInBatch(myPeakList, AnnotationData=plusAnno)
			TSSlength =temp$end_position - temp$start_position
			downstream = end(temp) - start(temp) + downstream
			temp$downstream=  downstream
			temp$TSSlength = TSSlength
	
			myList3 = as.data.frame(temp)
			rm(temp)
			
		if (dim(myList3)[1] != 0)
		{
			l3 =  cbind(as.character(myList3$feature), as.numeric(as.character(myList3$distancetoFeature)), as.numeric(rep(upstream, dim(myList3)[1])), as.numeric(as.character(myList3$downstream)), as.numeric(as.character(myList3$start_position)), as.numeric(as.character(myList3$end_position)))
			r3 = apply(l3, 1, getGeneSeq,genome)
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
		colnames(r)= c("feature", "distancetoFeature", "upstream", "downstream", "seq")
		r4 = merge(r, myList3)
		RangedData(IRanges(start=r4$start, end = r4$end, names =as.character(r4$name)), space=as.character(r4$space), upstream=rep(upstream, dim(r4)[1]), downstream= rep(downstream.bk, dim(r4)[1]), sequence=unlist(r4$seq))
	}
	else
	{
		stop("genome needs to be either a BSgenome object or Mart object!")
	}
}