findVennCounts <- function(Peaks, NameOfPeaks, maxgap=0L, minoverlap = 1L,  totalTest, useFeature=FALSE)
{
	if (missing(totalTest))
	{
		stop("Missing totalTest which is required! The parameter totalTest is the total number of possible Peaks in the testing space.")
	}
	if (missing(NameOfPeaks) ||  mode(NameOfPeaks) != "character")
	{
		stop("Missing required character vector NameOfPeaks")
	}
	if (missing(Peaks))
	{
		stop("Missing Peaks which is a list of peaks in RangedData!")
	}
	n1 = length(Peaks)
	n2 = length(NameOfPeaks)
	if (n1 <2)
	{
		stop("The number of element in NameOfPeaks is less than 1, need at least 2 peak ranges!")
	}
	if (n1 > n2)
	{
		stop("The number of element in NameOfPeaks is less than the number of elements in Peaks, need to be equal!")
	}
	if (n1 < n2)
	{
		stop("The number of element in NameOfPeaks is larger than the number of elements in Peaks, need to be equal!")
	}
	lapply(seq_len(n1), function(i)
	{
		if (class(Peaks[[i]]) != "RangedData")
		{
			err.msg = paste("Element", i, "in Peaks is not a valid RangedData object", sep=" ")
			stop(err.msg)
		}
	})
	lapply(seq_len(n1), function(i)
	{	
		if (!length(rownames(Peaks[i])))
		{
			rownames(Peaks[[i]]) = formatC(1:dim(Peaks[[i]])[1], width=nchar(dim(Peaks[[i]])[1]), flag='0')
		}		
	})
	if (n1 == 2)
	{
		a1 = cbind(c(0,0,1,1),c(0,1,0,1))
		colnames(a1) = NameOfPeaks
		a2 = vennCounts(a1)
		if (!useFeature)
		{
			overlappingPeaks = findOverlappingPeaks(Peaks[[1]],Peaks[[2]],NameOfPeaks1 = NameOfPeaks[1], NameOfPeaks2=NameOfPeaks[2], maxgap=maxgap,  minoverlap = minoverlap, select="first")
			if (length(overlappingPeaks$Peaks1withOverlaps) ==0)
			{
				p1.and.p2 = 0
			}
			else
			{	
				p1.inBoth = dim(overlappingPeaks$Peaks1withOverlaps)[1]
				p2.inBoth = dim(overlappingPeaks$Peaks2withOverlaps)[1]
				p1.and.p2 = min(p1.inBoth, p2.inBoth)
			}	
		#	p1.and.p2 = length(unique(overlappingPeaks$OverlappingPeaks[[NameOfPeaks[1]]]))
        #	print(overlappingPeaks)
			p1 = dim(Peaks[[1]])[1]
			p2 = dim(Peaks[[2]])[1]
        #	print(p1)
        #	print(p2)
			
		}
		else
		{
			if (length(Peaks[[1]]$feature) ==0 || length(Peaks[[2]]$feature) ==0)
			{
				warning("feature field in at least one of the Peaks RangedData has 0 elements!")				
			}
			p1.and.p2 = length(intersect(Peaks[[1]]$feature, Peaks[[2]]$feature))
			p1 = length(unique(Peaks[[1]]$feature))
			p2 = length(unique(Peaks[[2]]$feature))	
		}
		p1only = p1 - p1.and.p2
		p2only = p2 - p1.and.p2
		neither = totalTest - (p1 + p2 - p1.and.p2) 
		Counts =c(neither,p2only,p1only,p1.and.p2)
		if (p1only <0 || p2only <0)
                {
                        warning("negative counts generated when multiple peaks overlap with one peak!")
                }

		a2[,3] = Counts
		p.value = phyper(p1.and.p2 -1, p2, totalTest-p2, p1, lower.tail = FALSE,log.p = FALSE)
		list(p.value=p.value, vennCounts = a2)
	}
	else
	{
		stop("Larger than 2 list is not implemented yet!") 
	}
}
