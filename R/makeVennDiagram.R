makeVennDiagram <-
function(Peaks, NameOfPeaks, maxgap=0,  totalTest, cex=1.5, counts.col="red", useFeature=FALSE)
{
	if (missing(totalTest))
	{
		stop("Missing totalTest which is required! The parameter totalTest is the total number of possible Peaks in the testing space.")
	}
	if (missing(Peaks))
	{
		stop("Missing Peaks which is a list of peaks in RangedData!")
	}
	if (missing(NameOfPeaks) ||  mode(NameOfPeaks) != "character")
	{
		stop("Missing required character vector NameOfPeaks")
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
		x = findVennCounts(Peaks=Peaks,NameOfPeaks=NameOfPeaks,maxgap=maxgap,totalTest=totalTest,useFeature=useFeature)
		a2 = x$vennCounts
		p.value = x$p.value
		vennDiagram(a2,names = NameOfPeaks, cex=cex, counts.col = counts.col)
		list(p.value=p.value,vennCounts=a2)
	}
	else if (n1 == 3)
	{
		a1 = cbind(c(0, 0, 0,0,1, 1,1,1), c(0, 0,1,1, 0,0,1, 1),c(0,1,0,1,0,1,0,1))
		colnames(a1) = NameOfPeaks
		a2 = vennCounts(a1)
		x1 = findVennCounts(list(Peaks[[1]],Peaks[[2]]), NameOfPeaks=NameOfPeaks[1:2], maxgap=maxgap,
						totalTest=totalTest,useFeature=useFeature)
		p.value.1vs2 = x1$p.value
		counts1 = x1$vennCounts
		p1.and.p2 = counts1[4,3]
		
		x2 = findVennCounts(list(Peaks[[1]],Peaks[[3]]), NameOfPeaks=c(NameOfPeaks[1], NameOfPeaks[3]), maxgap=maxgap,
						totalTest=totalTest,useFeature=useFeature)
		p.value.1vs3 = x2$p.value
		counts2 = x2$vennCounts	
		p1.and.p3 = counts2[4,3]
		
		x3 = findVennCounts(list(Peaks[[2]],Peaks[[3]]), NameOfPeaks=NameOfPeaks[2:3], maxgap=maxgap,
						totalTest=totalTest,useFeature=useFeature)
		p.value.2vs3 = x3$p.value
		counts3 = x3$vennCounts
		p2.and.p3 = counts3[4,3]
		
		if(!useFeature)
		{
			overlappingPeaks123 = findOverlappingPeaks(
				findOverlappingPeaks(Peaks[[1]],Peaks[[2]],NameOfPeaks1 = NameOfPeaks[1], NameOfPeaks2=NameOfPeaks[2], maxgap=maxgap, select="first")$Peaks1withOverlap,
				Peaks[[3]], NameOfPeaks1 = NameOfPeaks[1], NameOfPeaks2=NameOfPeaks[3],maxgap=maxgap, select="first")$OverlappingPeaks
			p1.and.p2.and.p3 = length(unique(overlappingPeaks123[[NameOfPeaks[1]]]))
		
			p1 = length(rownames(Peaks[[1]]))
			p2 = length(rownames(Peaks[[2]]))
			p3 = length(rownames(Peaks[[3]]))
		}
		else
		{
			if (length(Peaks[[1]]$feature) ==0 || length(Peaks[[2]]$feature) ==0 ||  length(Peaks[[3]]$feature) ==0)
			{
				warning("feature field in at least one of the Peaks RangedData has 0 elements!")				
			}
			p1.and.p2.and.p3 = length(intersect(intersect(Peaks[[1]]$feature, Peaks[[2]]$feature),Peaks[[3]]$feature))
			p1 = length(unique(Peaks[[1]]$feature))
			p2 = length(unique(Peaks[[2]]$feature))	
			p3 = length(unique(Peaks[[3]]$feature))
		}
		p1only = p1 - p1.and.p2 - p1.and.p3 + p1.and.p2.and.p3
		p2only = p2 - p1.and.p2 - p2.and.p3 + p1.and.p2.and.p3
		p3only = p3 - p2.and.p3 - p1.and.p3 + p1.and.p2.and.p3
		if (p1only <0 || p2only <0 || p3only <0)
		{
			warning("negative counts generated when multiple peaks overlap with one peak!")
		}	
		neither123 = totalTest -p1only - p2only - p3only - p1.and.p2 - p1.and.p3 - p2.and.p3 + 2 * p1.and.p2.and.p3
		
		Counts =c(neither123,p3only,p2only,p2.and.p3-p1.and.p2.and.p3, p1only, p1.and.p3-p1.and.p2.and.p3, p1.and.p2-p1.and.p2.and.p3,p1.and.p2.and.p3)
		a2[,4] = Counts		
		vennDiagram(a2, names = NameOfPeaks)
		list(p.value.1vs2 = p.value.1vs2, p.value.1vs3 = p.value.1vs3, p.value.2vs3 = p.value.2vs3, vennCounts=a2)
	}
	else if(n1 == 4)
        {
		a1 = cbind(c(0, 0, 0,0,0, 0, 0,0,1, 1,1,1,1, 1,1,1), c(0, 0,1,1,0, 0,1,1,0, 0,1,1, 0,0,1, 1),c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1), c(0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1))
		colnames(a1) = NameOfPeaks
		a2 = vennCounts(a1)
		x1 = findVennCounts(list(Peaks[[1]],Peaks[[2]]), NameOfPeaks=NameOfPeaks[1:2], maxgap=maxgap,
						totalTest=totalTest,useFeature=useFeature)
		p.value.1vs2 = x1$p.value
		counts1 = x1$vennCounts
		p1.and.p2 = counts1[4,3]
		
		x2 = findVennCounts(list(Peaks[[1]],Peaks[[3]]), NameOfPeaks=c(NameOfPeaks[1], NameOfPeaks[3]), maxgap=maxgap,
						totalTest=totalTest,useFeature=useFeature)
		p.value.1vs3 = x2$p.value
		counts2 = x2$vennCounts	
		p1.and.p3 = counts2[4,3]
		
		x3 = findVennCounts(list(Peaks[[2]],Peaks[[3]]), NameOfPeaks=NameOfPeaks[2:3], maxgap=maxgap,
						totalTest=totalTest,useFeature=useFeature)
		p.value.2vs3 = x3$p.value
		counts3 = x3$vennCounts
		p2.and.p3 = counts3[4,3]
		
		x4 = findVennCounts(list(Peaks[[3]],Peaks[[4]]), NameOfPeaks=NameOfPeaks[3:4], maxgap=maxgap,
						totalTest=totalTest,useFeature=useFeature)
		p.value.3vs4 = x4$p.value
		counts4 = x4$vennCounts
		p3.and.p4 = counts4[4,3]
		
		x5 = findVennCounts(list(Peaks[[1]],Peaks[[4]]), NameOfPeaks=c(NameOfPeaks[1],NameOfPeaks[4]), maxgap=maxgap,totalTest=totalTest,useFeature=useFeature)
		p.value.1vs4 = x5$p.value
		counts5 = x5$vennCounts
		p1.and.p4 = counts5[4,3]
		
		x6 = findVennCounts(list(Peaks[[2]],Peaks[[4]]), NameOfPeaks=c(NameOfPeaks[2],NameOfPeaks[4]), maxgap=maxgap,totalTest=totalTest,useFeature=useFeature)
		p.value.2vs4 = x6$p.value
		counts6 = x6$vennCounts
		p2.and.p4 = counts6[4,3]		
		
		if (!useFeature)
		{
			overlappingPeaks123 = findOverlappingPeaks(
				findOverlappingPeaks(Peaks[[1]],Peaks[[2]],NameOfPeaks1 = NameOfPeaks[1], NameOfPeaks2=NameOfPeaks[2], maxgap=maxgap, select="first")$Peaks1withOverlap,
				Peaks[[3]], NameOfPeaks1 = NameOfPeaks[1], NameOfPeaks2=NameOfPeaks[3],maxgap=maxgap, select="first")
			p1.and.p2.and.p3 = length(unique(overlappingPeaks123$OverlappingPeaks[[NameOfPeaks[1]]]))
		
			overlappingPeaks124 = findOverlappingPeaks(
				findOverlappingPeaks(Peaks[[1]],Peaks[[2]],NameOfPeaks1 = NameOfPeaks[1], NameOfPeaks2=NameOfPeaks[2], maxgap=maxgap, select="first")$Peaks1withOverlap,
				Peaks[[4]], NameOfPeaks1 = NameOfPeaks[1], NameOfPeaks2=NameOfPeaks[4],maxgap=maxgap, select="first")$OverlappingPeaks
			p1.and.p2.and.p4 = length(unique(overlappingPeaks124[[NameOfPeaks[1]]]))
		
			overlappingPeaks234 = findOverlappingPeaks(
				findOverlappingPeaks(Peaks[[2]],Peaks[[3]],NameOfPeaks1 = NameOfPeaks[2], NameOfPeaks2=NameOfPeaks[3], maxgap=maxgap, select="first")$Peaks1withOverlap,
				Peaks[[4]], NameOfPeaks1 = NameOfPeaks[2], NameOfPeaks2=NameOfPeaks[4],maxgap=maxgap, select="first")$OverlappingPeaks
			p2.and.p3.and.p4 = length(unique(overlappingPeaks234[[NameOfPeaks[2]]]))
		
			overlappingPeaks134 = findOverlappingPeaks(
				findOverlappingPeaks(Peaks[[1]],Peaks[[3]],NameOfPeaks1 = NameOfPeaks[1], NameOfPeaks2=NameOfPeaks[3], maxgap=maxgap, select="first")$Peaks1withOverlap,
				Peaks[[4]], NameOfPeaks1 = NameOfPeaks[1], NameOfPeaks2=NameOfPeaks[4],maxgap=maxgap, select="first")$OverlappingPeaks
			p1.and.p3.and.p4 = length(unique(overlappingPeaks134[[NameOfPeaks[1]]]))
				
			overlappingPeaks1234 = findOverlappingPeaks(overlappingPeaks123$Peaks1withOverlap,
				Peaks[[4]], NameOfPeaks1 = NameOfPeaks[1], NameOfPeaks2=NameOfPeaks[4],maxgap=maxgap, select="first")$OverlappingPeaks
			p1.and.p2.and.p3.and.p4 = length(unique(overlappingPeaks1234[[NameOfPeaks[1]]]))
		
			p1 = length(rownames(Peaks[[1]]))
			p2 = length(rownames(Peaks[[2]]))
			p3 = length(rownames(Peaks[[3]]))
			p4 = length(rownames(Peaks[[4]]))
		}
		else
		{
			if (length(Peaks[[1]]$feature) ==0 || length(Peaks[[2]]$feature) ==0 ||  length(Peaks[[3]]$feature) ==0 || length(Peaks[[4]]$feature) ==0)
			{
				warning("feature field in at least one of the Peaks RangedData has 0 elements!")				
			}
			p1.and.p2.and.p3 = length(intersect(intersect(Peaks[[1]]$feature, Peaks[[2]]$feature),Peaks[[3]]$feature))
			p1.and.p2.and.p4 = length(intersect(intersect(Peaks[[1]]$feature, Peaks[[2]]$feature),Peaks[[4]]$feature))
			p2.and.p3.and.p4 = length(intersect(intersect(Peaks[[2]]$feature, Peaks[[3]]$feature),Peaks[[4]]$feature))
			p1.and.p3.and.p4 = length(intersect(intersect(Peaks[[1]]$feature, Peaks[[3]]$feature),Peaks[[4]]$feature))
			p1.and.p2.and.p3.and.p4 = length(intersect(intersect(intersect(Peaks[[1]]$feature, Peaks[[2]]$feature),Peaks[[3]]$feature),Peaks[[4]]$feature))
			
			p1 = length(unique(Peaks[[1]]$feature))
			p2 = length(unique(Peaks[[2]]$feature))	
			p3 = length(unique(Peaks[[3]]$feature))
			p4 = length(unique(Peaks[[4]]$feature))
		}

		p1only = p1 - p1.and.p2 - p1.and.p3 + p1.and.p2.and.p3 - p1.and.p4 + p1.and.p2.and.p4 + p1.and.p3.and.p4 - p1.and.p2.and.p3.and.p4
		p2only =  p2 - p1.and.p2 - p2.and.p3 + p1.and.p2.and.p3 - p2.and.p4 + p2.and.p3.and.p4 + p1.and.p2.and.p4 - p1.and.p2.and.p3.and.p4
		p3only =  p3 - p1.and.p3 - p2.and.p3 + p1.and.p2.and.p3 - p3.and.p4 + p2.and.p3.and.p4 + p1.and.p3.and.p4 - p1.and.p2.and.p3.and.p4
		p4only =  p4 - p1.and.p4 - p2.and.p4 + p1.and.p2.and.p4 - p3.and.p4 + p2.and.p3.and.p4 + p1.and.p3.and.p4 - p1.and.p2.and.p3.and.p4
	
		 if (p1only <0 || p2only <0 || p3only <0 || p4only <0)
                {
                        warning("negative counts generated when multiple peaks overlap with one peak!")
                }
	
		neither1234 = totalTest - ( p1 + p2 - p1.and.p2 + p3only + p4only + p3.and.p4
				- p1.and.p3.and.p4 - p2.and.p3.and.p4 + p1.and.p2.and.p3.and.p4)
		Counts =c(neither1234, p4only, p3only, p3.and.p4 - p1.and.p3.and.p4 - p2.and.p3.and.p4 + p1.and.p2.and.p3.and.p4, p2only, p2.and.p4 - p2.and.p3.and.p4 - p1.and.p2.and.p4 + p1.and.p2.and.p3.and.p4, p2.and.p3 - p2.and.p3.and.p4 - p1.and.p2.and.p3 + p1.and.p2.and.p3.and.p4, 
p2.and.p3.and.p4 - p1.and.p2.and.p3.and.p4, p1only,
		p1.and.p4 - p1.and.p3.and.p4 - p1.and.p2.and.p4 + p1.and.p2.and.p3.and.p4,
		p1.and.p3 - p1.and.p3.and.p4 - p1.and.p2.and.p3 + p1.and.p2.and.p3.and.p4,
		p1.and.p3.and.p4 - p1.and.p2.and.p3.and.p4,
		p1.and.p2 - p1.and.p2.and.p3 - p1.and.p2.and.p4 + p1.and.p2.and.p3.and.p4,
		p1.and.p2.and.p4 - p1.and.p2.and.p3.and.p4,
		p1.and.p2.and.p3 - p1.and.p2.and.p3.and.p4,
		p1.and.p2.and.p3.and.p4
		)
		a2[,5] = Counts	
		a = cbind(a2[,5], a2[,1:4])
		colnames(a)[1]="num"
		#vennDiagram(a2, names = NameOfPeaks)
		rownames(a)=c("0000", "0001", "0010", "0011", "0100", "0101", "0110", "0111", "1000", "1001", "1010",
		 "1011", "1100", "1101", "1110", "1111")
		gplots:::drawVennDiagram(a)
		list(p.value.1vs2 = p.value.1vs2, p.value.1vs3 = p.value.1vs3, p.value.2vs3 = p.value.2vs3, p.value.1vs4= p.value.1vs4, p.value.2vs4= p.value.2vs4, p.value.3vs4= p.value.3vs4, vennCounts=a2)
	}	
	else
	{
		stop("Larger than 4 lists are not implemented yet")
	}
}
