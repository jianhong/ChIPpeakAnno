findOverlappingPeaks <-
function(Peaks1,Peaks2,maxgap=100,multiple=c(TRUE, FALSE), NameOfPeaks1="TF1", NameOfPeaks2="TF2")
{
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
		if (!length(rownames(Peaks1)))
		{
			rownames(Peaks1) = formatC(1:dim(Peaks1)[1], width=nchar(dim(Peaks1)[1]), flag='0')
		}		
		if (!length(rownames(Peaks2)))
		{
			rownames(Peaks2) = formatC(1:dim(Peaks2)[1], width=nchar(dim(Peaks2)[1]), flag='0')
		}

		r2 = cbind(rownames(Peaks2), start(Peaks2), end(Peaks2),Peaks2$strand)
		colnames(r2) = c(NameOfPeaks2, paste(NameOfPeaks2, "start", sep="_"), paste(NameOfPeaks2,"end", sep="_"),"strand")
		allChr.Anno = unique(space(Peaks2))
			
		numberOfChromosome = length(unique(space(Peaks1)))
		allChr = unique(as.character(space(Peaks1)))
		allChr = sub(' +', '',allChr)
		
		if(length(grep("chr", allChr, fixed=TRUE))>0 && length(grep("chr", allChr.Anno, fixed=TRUE))==0)
		{
			allChr = sub("chr", '', allChr)
			Peaks1 = RangedData(IRanges(start=start(Peaks1), 
							end= end(Peaks1),
							names = rownames(Peaks1)),
							space = sub("chr", "", space(Peaks1)),
							strand = Peaks1$strand
							)
		}
		
		if(length(grep("chr", allChr, fixed=TRUE))==0 && length(grep("chr", allChr.Anno, fixed=TRUE))>0)
		{
			allChr = paste("chr", allChr, sep="")
			Peaks1 = RangedData(IRanges(start=start(Peaks1), 
							end= end(Peaks1),
							names = rownames(Peaks1)),
							space = paste("chr", space(Peaks1), sep=""),
							strand = Peaks1$strand
							)
		}
		
		z1 = cbind(as.character(rownames(Peaks1)), as.character(space(Peaks1)),start(Peaks1), end(Peaks1), Peaks1$strand)
		if (length(Peaks1$strand) ==0)
		{
			z1 = cbind(z1,rep("1", dim(z1)[1]))
		}
		colnames(z1) = c(NameOfPeaks1, "chr", paste(NameOfPeaks1, "start", sep="_"), paste(NameOfPeaks1,"end", sep="_"),"strand1")
		z1[,2] = sub(' +', '',z1[,2])
		
	   r1 = do.call(rbind, lapply(seq_len(numberOfChromosome), function(i) 	{
       chr = allChr[i]
		  	if(chr %in% allChr.Anno)
			{
				Ranges1 = sort(IRanges(start=start(Peaks1[chr]),end=end(Peaks1[chr]), names=rownames(Peaks1[chr])))
				Ranges2 = sort(IRanges(start=start(Peaks2[chr]),end=end(Peaks2[chr]), names=rownames(Peaks2[chr])))
				tree = IntervalTree(Ranges2)
				matches = findOverlaps(tree, query = Ranges1, maxgap = maxgap, multiple=multiple)
				if (multiple)
				{
					matchmatrix = matchMatrix(matches)
					qname = names(Ranges1)[matchmatrix[,1]]
					tname = names(Ranges2)[matchmatrix[,2]]
					data.frame(Peaks1 = qname,
                    		chr = rep(chr,dim(matchmatrix)[1]),
                    		Peaks2 =  tname)
				}
				else
				{
					qname = names(Ranges1)
					tname = names(Ranges2)[matches]
					data.frame(Peaks1 = qname,
                    		chr = rep(chr,length(qname)),
                    		Peaks2 =  tname)
				}
			}
       }))
		colnames(r1) = c(NameOfPeaks1, "chr", NameOfPeaks2) 
		r3 = merge(r1,r2, by=NameOfPeaks2)		
		r = merge(r3, z1)
		shortestDistance = apply(cbind(abs(as.numeric(as.character(r[[paste(NameOfPeaks2, "start", sep="_")]]))-as.numeric(as.character(r[[paste(NameOfPeaks1, "end", sep="_")]]))), 
			abs(as.numeric(as.character(r[[paste(NameOfPeaks2, "start", sep="_")]]))-as.numeric(as.character(r[[paste(NameOfPeaks1, "start", sep="_")]]))), 
			abs(as.numeric(as.character(r[[paste(NameOfPeaks2, "end", sep="_")]]))-as.numeric(as.character(r[[paste(NameOfPeaks1, "end", sep="_")]]))),
			abs(as.numeric(as.character(r[[paste(NameOfPeaks2, "end", sep="_")]]))-as.numeric(as.character(r[[paste(NameOfPeaks1, "start", sep="_")]]))))
			,1,min)	
		overlapFeature = do.call(c, lapply(seq(from=1,to=dim(r)[1],by=1), function(i) 	{
		if (as.numeric(as.character(r[[paste(NameOfPeaks1, "start", sep="_")]][i]))>= as.numeric(as.character(r[[paste(NameOfPeaks2, "start", sep="_")]][i])) &
			as.numeric(as.character(r[[paste(NameOfPeaks1, "start", sep="_")]][i]))<= as.numeric(as.character(r[[paste(NameOfPeaks2, "end", sep="_")]][i])))
		{
			if (as.numeric(as.character(r[[paste(NameOfPeaks1, "end", sep="_")]][i])) >= as.numeric(as.character(r[[paste(NameOfPeaks2, "start", sep="_")]][i])) &
			as.numeric(as.character(r[[paste(NameOfPeaks1, "end", sep="_")]][i])) <= as.numeric(as.character(r[[paste(NameOfPeaks2, "end", sep="_")]][i])))
			{
				"inside"
			}
			else if (as.character(r$strand[i]) == "1" |  as.character(r$strand[i]) =="+")
			{
				"overlapEnd"
			}
			else
			{
				"overlapStart"
			}
		}
		else if (as.numeric(as.character(r[[paste(NameOfPeaks1, "end", sep="_")]][i])) >= as.numeric(as.character(r[[paste(NameOfPeaks2, "start", sep="_")]][i])) &
			as.numeric(as.character(r[[paste(NameOfPeaks1, "end", sep="_")]][i])) <= as.numeric(as.character(r[[paste(NameOfPeaks2, "end", sep="_")]][i])))
		{
			if (as.character(r$strand[i]) == "1" |  as.character(r$strand[i]) =="+")
			{
				"overlapStart"
			}
			else
			{
				"overlapEnd"
			}
		}
		else if (as.numeric(as.character(r[[paste(NameOfPeaks1, "end", sep="_")]][i])) >= as.numeric(as.character(r[[paste(NameOfPeaks2, "end", sep="_")]][i])) &
			as.numeric(as.character(r[[paste(NameOfPeaks1, "start", sep="_")]][i])) <= as.numeric(as.character(r[[paste(NameOfPeaks2, "start", sep="_")]][i])))
		{
				 "includeFeature"
		}
		else if (as.numeric(as.character(r[[paste(NameOfPeaks1, "end", sep="_")]][i])) < as.numeric(as.character(r[[paste(NameOfPeaks2, "start", sep="_")]][i])))
		{
			if (as.character(r$strand[i]) == "1" |  as.character(r$strand[i]) =="+")
			{
				"upstream"
			}
			else
			{
				"downstream"
			}
		}
		else if (as.numeric(as.character(r[[paste(NameOfPeaks1, "start", sep="_")]][i])) > as.numeric(as.character(r[[paste(NameOfPeaks2, "end", sep="_")]][i])))
		{
			if (as.character(r$strand[i]) == "1" |  as.character(r$strand[i]) =="+")
			{
				"downstream"
			}
			else
			{
				"upstream"
			}
		}
		}))
		
		r$overlapFeature = overlapFeature
		r$shortestDistance = shortestDistance
		minstart = apply(cbind(as.numeric(as.character(r[[paste(NameOfPeaks1,"start", sep="_")]])), as.numeric(as.character(r[[paste(NameOfPeaks2,"start", sep="_")]]))), MARGIN=1, min)
		maxend = apply(cbind(as.numeric(as.character(r[[paste(NameOfPeaks1,"end", sep="_")]])), as.numeric(as.character(r[[paste(NameOfPeaks2,"end", sep="_")]]))), MARGIN =1, max)
		MergedPeaks = RangedData(IRanges(start = minstart, end = maxend,
			names = paste(NameOfPeaks1,as.character(r[[NameOfPeaks1]]), NameOfPeaks2, as.character(r[[NameOfPeaks2]]), sep="-")), 
			space=as.character(r$chr))	
		#MergedPeaks <- MergedPeaks[order(rownames(MergedPeaks)),]
		r1 = unique(cbind(as.character(r[[NameOfPeaks1]]),
					as.character(r$chr),
					as.numeric(as.character(r[[paste(NameOfPeaks1,"start",sep="_")]])),
					as.numeric(as.character(r[[paste(NameOfPeaks1,"end",sep="_")]])),
					as.character(r$strand1))
					)
		r2 = unique(cbind(as.character(r[[NameOfPeaks2]]),
					as.character(r$chr),
					as.numeric(as.character(r[[paste(NameOfPeaks2,"start",sep="_")]])),
					as.numeric(as.character(r[[paste(NameOfPeaks2,"end",sep="_")]])),
					as.character(r$strand))
					)
		Peaks1withOverlaps = RangedData(IRanges(start=as.numeric(as.character(r1[,3])), 
							end= as.numeric(as.character(r1[,4])),
							names = as.character(r1[,1])),
							space = as.character(r1[,2]), 
							strand=as.character(r1[,5]))
		Peaks2withOverlaps = RangedData(IRanges(start=as.numeric(as.character(r2[,3])), 
							end= as.numeric(as.character(r2[,4])),
							names = as.character(r2[,1])),
							space = as.character(r2[,2]), 
							strand=as.character(r2[,5]))
		list(OverlappingPeaks = r[order(r[[NameOfPeaks1]]),], MergedPeaks = MergedPeaks,Peaks1withOverlaps=Peaks1withOverlaps, Peaks2withOverlaps=Peaks2withOverlaps)
	}

