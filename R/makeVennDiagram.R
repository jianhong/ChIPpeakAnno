makeVennDiagram <-
function(Peaks1, Peaks2, overlappingPeaks,totalTest, NameOfPeaks1="TF1", NameOfPeaks2="TF2", cex=1, counts.col="red")
{
	if (missing(totalTest))
	{
		stop("Missing totalTest which is required! The parametr totalTest is the total number of possible Peaks in the testing space.")
	}
	if (missing(Peaks1))
		{
			stop("Missing Peaks1 which is required!")
		}
		if (class(Peaks1) != "RangedData")
		{
			stop("No valid Peaks1 passed in. It needs to be RangedData object")
		}
		if (missing(Peaks2))
		{
			stop("No Peaks2 as RangedData is passed in")
		}
		if (class(Peaks2) != "RangedData")
		{
			stop("Peaks2 needs to be RangedData object")
		}	
		if (missing(overlappingPeaks))
		{
			stop("No overlappingPeaks as dataframe is passed in, which can be obtained with FindOverlappingPeaks(Peaks1,Peaks2, multiple=F)$OverlappingPeaks")
		}	
		if (!length(rownames(Peaks1)))
		{
			rownames(Peaks1) = formatC(1:dim(Peaks1)[1], width=nchar(dim(Peaks1)[1]), flag='0')
		}		
		if (!length(rownames(Peaks2)))
		{
			rownames(Peaks2) = formatC(1:dim(Peaks2)[1], width=nchar(dim(Peaks2)[1]), flag='0')
		}

	a1 = cbind(c(0,0,1,1),c(0,1,0,1))
	colnames(a1) = c(NameOfPeaks1, NameOfPeaks2)
	a2 = vennCounts(a1)
	p1.and.p2 = length(unique(overlappingPeaks[[NameOfPeaks1]]))
	p1 = length(rownames(Peaks1))
	p2 = length(rownames(Peaks2))
	p1only = length(setdiff(as.character(rownames(Peaks1)), as.character(overlappingPeaks[[NameOfPeaks1]])))
	p2only = length(setdiff(as.character(rownames(Peaks2)), as.character(overlappingPeaks[[NameOfPeaks2]])))
	p2inp1 = length(unique(overlappingPeaks[[NameOfPeaks2]]))
	neither = totalTest - p1only - p2 
	Counts =c(neither,p2only,p1only,p1.and.p2)
	a2[,3] = Counts
	vennDiagram(a2,names = c(NameOfPeaks1, NameOfPeaks2), cex=cex, counts.col = counts.col)
	p.value = phyper(p1.and.p2 -1, p2, totalTest-p2, p1, lower.tail = FALSE,log.p = FALSE)
	list(p.value=p.value)
}

