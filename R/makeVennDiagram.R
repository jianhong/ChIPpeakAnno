makeVennDiagram <-
function(Peaks, NameOfPeaks, maxgap=0,  totalTest, cex=1.5, counts.col="red")
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
		x = findVennCounts(Peaks=Peaks,NameOfPeaks=NameOfPeaks,maxgap=maxgap,totalTest=totalTest)
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
						totalTest=totalTest)
		p.value.1vs2 = x1$p.value
		counts1 = x1$vennCounts
		p1.and.p2 = counts1[4,3]
		
		x2 = findVennCounts(list(Peaks[[1]],Peaks[[3]]), NameOfPeaks=c(NameOfPeaks[1], NameOfPeaks[3]), maxgap=maxgap,
						totalTest=totalTest)
		p.value.1vs3 = x2$p.value
		counts2 = x2$vennCounts	
		p1.and.p3 = counts2[4,3]
		
		x3 = findVennCounts(list(Peaks[[2]],Peaks[[3]]), NameOfPeaks=NameOfPeaks[2:3], maxgap=maxgap,
						totalTest=totalTest)
		p.value.2vs3 = x3$p.value
		counts3 = x3$vennCounts
		p2.and.p3 = counts3[4,3]
		
		overlappingPeaks123 = findOverlappingPeaks(
				findOverlappingPeaks(Peaks[[1]],Peaks[[2]],NameOfPeaks1 = NameOfPeaks[1], NameOfPeaks2=NameOfPeaks[2], maxgap=maxgap, multiple=F)$Peaks1withOverlap,
				Peaks[[3]], NameOfPeaks1 = NameOfPeaks[1], NameOfPeaks2=NameOfPeaks[3],maxgap=maxgap, multiple=F)$OverlappingPeaks
		p1.and.p2.and.p3 = length(unique(overlappingPeaks123[[NameOfPeaks[1]]]))
		
		p1 = length(rownames(Peaks[[1]]))
		p2 = length(rownames(Peaks[[2]]))
		p3 = length(rownames(Peaks[[3]]))
		p1only = p1 - p1.and.p2 - p1.and.p3 + p1.and.p2.and.p3
		p2only = p2 - p1.and.p2 - p2.and.p3 + p1.and.p2.and.p3
		p3only = p3 - p2.and.p3 - p1.and.p3 + p1.and.p2.and.p3
		
		neither123 = totalTest -p1only - p2only - p3only - p1.and.p2 - p1.and.p3 - p2.and.p3 + 2 * p1.and.p2.and.p3
		
		Counts =c(neither123,p3only,p2only,p2.and.p3, p1only, p1.and.p3, p1.and.p2,p1.and.p2.and.p3)
		a2[,4] = Counts		
		vennDiagram(a2, names = NameOfPeaks)
		p.value.overall3 = phyper(p1.and.p2.and.p3 -1, p3, totalTest-p3, p1.and.p2, lower.tail = FALSE,log.p = FALSE)	
		p.value.overall1 = phyper(p1.and.p2.and.p3 -1, p1, totalTest-p1, p2.and.p3, lower.tail = FALSE,log.p = FALSE)
  		p.value.overall2 = phyper(p1.and.p2.and.p3 -1, p2, totalTest-p2, p1.and.p3, lower.tail = FALSE,log.p = FALSE)
		p.value.overall = max(p.value.overall1, p.value.overall2, p.value.overall3)
		list(p.value.1vs2 = p.value.1vs2, p.value.1vs3 = p.value.1vs3, p.value.2vs3 = p.value.2vs3, vennCounts=a2)
	}
	else
	{
		stop("Larger than 3 lists are not implemented yet")
	}
}

#source("~/dev/svndev/ChIPpeakAnno/R/makeVennDiagram.R")
#peaks1 = RangedData(IRanges(start = c(967654, 2010897, 2496704), end = c(967754, 2010997, 2496804), names = c("Site1", "Site2", "Site3")), space = c("1", "2", "3"), strand=as.integer(1))
#peaks2 = RangedData(IRanges(start = c(967659, 2010898,  2496700, 3075866, 3123260), end = c(967869, 2011108, 2496920, 3076166, 3123470), names = c("t1", "t2", "t3", "t4", "t5")), space = c("1", "2", "3", "1", "2"), strand = c(1, 1, -1,-1,1))
#makeVennDiagram(list(peaks1,peaks2,peaks2),NameOfPeaks,totalTest=100)
#makeVennDiagram(list(peaks1,peaks1,peaks1),NameOfPeaks,totalTest=100)
#makeVennDiagram(list(peaks1,peaks2,peaks1),NameOfPeaks,totalTest=100)
