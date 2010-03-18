annotatePeakInBatch <-
function(myPeakList, mart,featureType=c("TSS","miRNA", "Exon"), AnnotationData,output=c("nearestStart", "overlapping","both"),multiple=c(FALSE,TRUE), maxgap=0)
{
		featureType = match.arg(featureType)
		if (missing(output))
		{
			output = "nearestStart"
		}
		if ((output == "overlapping" || output == "both") && missing(multiple))
		{
			stop("Missing requried logical argument multiple. It is requried when output is overlapping or both!") 
		}
		
		if (missing(myPeakList))
		{
			stop("Missing required argument myPeakList!")
		}
		if (class(myPeakList) != "RangedData")
		{
			stop("No valid myPeakList passed in. It needs to be RangedData object")
		}
		if (missing(AnnotationData))
		{
			message("No AnnotationData as RangedData is passed in, so now querying biomart database for AnnotationData ....")
			if (missing(mart) || class(mart) !="Mart")
			{
				stop("Error in querying biomart database. No valid mart object is passed in! Suggest call getAnnotation before calling annotatePeakInBatch")
			}
			AnnotationData<- getAnnotation(mart, feature=featureType)
			message("Done querying biomart database, start annotating ....Better way would be calling getAnnotation before calling annotatePeakInBatch")
		}
		if (class(AnnotationData) != "RangedData")
		{
			stop("AnnotationData needs to be RangedData object")
		}
		
		TSS.ordered <-AnnotationData
		rm(AnnotationData)
		if (!length(rownames(TSS.ordered)))
		{
			rownames(TSS.ordered) = formatC(1:dim(TSS.ordered)[1], width=nchar(dim(TSS.ordered)[1]), flag='0')
		}
		r2 = cbind(rownames(TSS.ordered), start(TSS.ordered), end(TSS.ordered), TSS.ordered$strand)
		colnames(r2) = c("feature_id", "start_position", "end_position", "strand")
		allChr.Anno = unique(space(TSS.ordered))
		
		#allChr.Anno = sub("chr", '', allChr.Anno)
		
		numberOfChromosome = length(unique(space(myPeakList)))
	
		if (!length(rownames(myPeakList)))
		{
			rownames(myPeakList) = formatC(1:dim(myPeakList)[1], width=nchar(dim(myPeakList)[1]), flag='0')
		}
		
		allChr = unique(as.character(space(myPeakList)))
		allChr = sub(' +', '',allChr)
		
		if(length(grep("chr", allChr, fixed=TRUE))>0 && length(grep("chr", allChr.Anno, fixed=TRUE))==0)
		{
			allChr = sub("chr", '', allChr)
			myPeakList = RangedData(IRanges(start=start(myPeakList), 
							end= end(myPeakList),
							names = rownames(myPeakList)),
							space = sub("chr", "", space(myPeakList))
							)
		}
		
		if(length(grep("chr", allChr, fixed=TRUE))==0 && length(grep("chr", allChr.Anno, fixed=TRUE))>0)
		{
			allChr = paste("chr", allChr, sep="")
			myPeakList = RangedData(IRanges(start=start(myPeakList), 
							end= end(myPeakList),
							names = rownames(myPeakList)),
							space = paste("chr", space(myPeakList), sep="")
							)
		}
		if (output == "nearestStart" || output == "both" || output == "n" || output == "b")
		{
			z1 = cbind(as.character(rownames(myPeakList)), as.character(space(myPeakList)),start(myPeakList), end(myPeakList))
			colnames(z1) = c("name", "chr", "peakStart", "peakEnd")
			z1[,2] = sub(' +', '',z1[,2])
		
			plusAnno = TSS.ordered[TSS.ordered$strand==1 | TSS.ordered$strand =="+",]
			minusAnno = TSS.ordered[TSS.ordered$strand== -1 | TSS.ordered$strand =="-",]
			r1 = do.call(rbind, lapply(seq_len(numberOfChromosome), function(i) {
				chr = allChr[i]
				if(chr %in% allChr.Anno)
				{
           		featureStart = c(start(plusAnno[chr]), end(minusAnno[chr]))				
				peakStart = start(myPeakList[chr])
				peakEnd = end(myPeakList[chr])
				name = rownames(myPeakList[chr])
				peakRanges = IRanges(start=peakStart, end=peakEnd, names=name)
				featureID  = c(rownames(plusAnno[chr]), rownames(minusAnno[chr]))
				featureRanges = IRanges(start=featureStart, end=featureStart, names=featureID)
				nearestFeature = featureRanges[nearest(peakRanges,featureRanges)]
         		data.frame(name = name,
                    chr = rep(chr,length(peakStart)),
                    peakStart = peakStart,
                    peakEnd = peakEnd,
                    feature_id =  names(nearestFeature))
				}
			}))
		
		if (length(r1)>0)
		{
			r3 = merge(r1,r2, by="feature_id")
			r = merge(r3, z1, all.y=TRUE)
			r.n = r
			r11 = r[!is.na(r$strand) & (r$strand==1 | r$strand =="+"),]
			r22 = r[!is.na(r$strand) & (r$strand==-1 | r$strand =="-"),]
			r33 = r[is.na(r$strand),]
			r33$insideFeature = replicate(length(r33$name), NA)
			r33$distancetoFeature = replicate(length(r33$name), NA)
		
		if (dim(r11)[1] >0)
		{
		distancetoFeature =  as.numeric(as.character(r11$peakStart)) - as.numeric(as.character(r11$start_position))
		length =  as.numeric(as.character(r11$end_position)) -  as.numeric(as.character(r11$start_position))
		#insideFeature = (distancetoFeature>=0 & distancetoFeature <= length)
		insideFeature = do.call(c, lapply(seq(from=1,to=dim(r11)[1],by=1), function(i) 	{
		if (as.numeric(as.character(r11$peakStart[i]))>= as.numeric(as.character(r11$start_position[i])) &
			as.numeric(as.character(r11$peakStart[i]))<= as.numeric(as.character(r11$end_position[i])))
		{
			if (as.numeric(as.character(r11$peakEnd[i])) >= as.numeric(as.character(r11$start_position[i])) &
			as.numeric(as.character(r11$peakEnd[i])) <= as.numeric(as.character(r11$end_position[i])))
			{
				"inside"
			}
			else
			{
				"overlapEnd"
			}
		}
		else if (as.numeric(as.character(r11$peakEnd[i])) >= as.numeric(as.character(r11$start_position[i])) &
			as.numeric(as.character(r11$peakEnd[i])) <= as.numeric(as.character(r11$end_position[i])))
		{
				"overlapStart"
		}
		else if (as.numeric(as.character(r11$peakEnd[i])) >= as.numeric(as.character(r11$end_position[i])) &
			as.numeric(as.character(r11$peakStart[i])) <= as.numeric(as.character(r11$start_position[i])))
		{
				 "includeFeature"
		}
		else if (as.numeric(as.character(r11$peakEnd[i])) < as.numeric(as.character(r11$start_position[i])))
		{
				"upstream"
		}
		else if (as.numeric(as.character(r11$peakStart[i])) > as.numeric(as.character(r11$end_position[i])))
		{
				"downstream"
		}
		}))
		r11$insideFeature = insideFeature		
		r11$distancetoFeature = distancetoFeature
		r.n= r11
		}
		if (dim(r22)[1] >0)
		{		
		distancetoFeature =  as.numeric(as.character(r22$end_position)) - as.numeric(as.character(r22$peakStart))
		length =  as.numeric(as.character(r22$end_position)) -  as.numeric(as.character(r22$start_position))
		#insideFeature = (distancetoFeature>=0 & distancetoFeature <= length)	
		insideFeature = do.call(c, lapply(seq(from=1,to=dim(r22)[1],by=1), function(i) 	{		
			if (as.numeric(as.character(r22$peakStart[i]))>= as.numeric(as.character(r22$start_position[i])) &
			as.numeric(as.character(r22$peakStart[i]))<= as.numeric(as.character(r22$end_position[i])))
			{
			if (as.numeric(as.character(r22$peakEnd[i])) >= as.numeric(as.character(r22$start_position[i])) &
			as.numeric(as.character(r22$peakEnd[i])) <= as.numeric(as.character(r22$end_position[i])))
			{
				"inside"
			}
			else
			{
				"overlapStart"
			}
			}
		else if (as.numeric(as.character(r22$peakEnd[i])) >= as.numeric(as.character(r22$start_position[i])) &
			as.numeric(as.character(r22$peakEnd[i])) <= as.numeric(as.character(r22$end_position[i])))
		{
				"overlapEnd"
		}
		else if (as.numeric(as.character(r22$peakEnd[i])) >= as.numeric(as.character(r22$end_position[i])) &
			as.numeric(as.character(r22$peakStart[i])) <= as.numeric(as.character(r22$start_position[i])))
		{
				 "includeFeature"
		}
		else if (as.numeric(as.character(r22$peakEnd[i])) < as.numeric(as.character(r22$start_position[i])))
		{
				"downstream"
		}
		else if (as.numeric(as.character(r22$peakStart[i])) > as.numeric(as.character(r22$end_position[i])))
		{
				"upstream"
		}
		}))					
							
		r22$insideFeature = insideFeature
		r22$distancetoFeature = distancetoFeature
		if (dim(r11)[1] >0)
		{
			r.n = rbind(r.n,r22)
		}
		else
		{
			r.n = r22
		}
		if (dim(r33)[1] >0)
		{
			if (dim(r11)[1] > 0 || dim(r22)[1] >0)
			{
				r.n = rbind(r.n, r33)
			}
			else
			{
				r.n = r33
			}
		}
		r.n$fromOverlappingOrNearest = rep("NearestStart", dim(r.n)[1])
		shortestDistance = apply(cbind(abs(as.numeric(as.character(r.n$start_position))-as.numeric(as.character(r.n$peakEnd))), 
					abs(as.numeric(as.character(r.n$start_position))-as.numeric(as.character(r.n$peakStart))), 
					abs(as.numeric(as.character(r.n$end_position))-as.numeric(as.character(r.n$peakEnd))),
					abs(as.numeric(as.character(r.n$end_position))-as.numeric(as.character(r.n$peakStart))))
		,1,min)	
		r.n$shortestDistance = shortestDistance
		
		} ## if dim(r22)[1] >0
		} ## if length(r1) >0
		} ## if output == "nearestStart"
		if (output == "overlapping" || output == "both" || output == "o" || output == "b")
		{
			r.o = findOverlappingPeaks(myPeakList,TSS.ordered ,maxgap=maxgap,multiple=multiple, NameOfPeaks1="peak", NameOfPeaks2="feature")$OverlappingPeaks
			r.o$fromOverlappingOrNearest = rep("Overlapping", dim(r.o)[1])
			distancetoFeature =  do.call(c, lapply(seq_len(dim(r.o)[1]), function(i)
			{
				if (as.character(r.o$strand[i])=="1" | as.character(r.o$strand[i]) =="+")
				{
					as.numeric(as.character(r.o$peak_start[i])) - as.numeric(as.character(r.o$feature_start[i]))
				}
				else
				{
					as.numeric(as.character(r.o$feature_end[i])) - as.numeric(as.character(r.o$peak_start[i]))
				}
			}))
			r.o$distancetoFeature = distancetoFeature
		}
		if (output == "nearestStart" || output == "n")
		{
			if(length(r1)>0)
			{
				r.output = RangedData(IRanges(start=as.numeric(as.character(r.n$peakStart)), 
						end=as.numeric(as.character(r.n$peakEnd)),
						names=paste(as.character(r.n$name),as.character(r.n$feature_id))),
          		peak = as.character(r.n$name),
				strand = as.character(r.n$strand), 
				feature = as.character(r.n$feature_id),
				start_position= as.numeric(as.character(r.n$start_position)),
				end_position=as.numeric(as.character(r.n$end_position)),
				insideFeature=as.character(r.n$insideFeature), 
				distancetoFeature=as.numeric(as.character(r.n$distancetoFeature)),
				shortestDistance = as.numeric(as.character(r.n$shortestDistance)),
				fromOverlappingOrNearest = as.character(r.n$fromOverlappingOrNearest),
				space = as.character(r.n$chr)
				)
				#r.output[order(rownames(r.output)),]
			}
			else
			{
				r.output = myPeakList
			}
		}
		else if ((output =="overlapping" || output =="o") || ((output =="both" || output=="b") && length(r1) ==0))
		{
			r.output=RangedData(IRanges(start=as.numeric(as.character(r.o$peak_start)), end=as.numeric(as.character(r.o$peak_end)),
								names=paste(as.character(r.o$peak),as.character(r.o$feature))), 
								peak = as.character(r.o$peak),
								strand = as.character(r.o$strand),
								feature = as.character(r.o$feature),
								start_position= as.numeric(as.character(r.o$feature_start)),
								end_position= as.numeric(as.character(r.o$feature_end)),
								insideFeature=r.o$overlapFeature,
								distancetoFeature=as.numeric(as.character(r.o$distancetoFeature)),
								shortestDistance = as.numeric(as.character(r.o$shortestDistance)),
								fromOverlappingOrNearest = as.character(r.o$fromOverlappingOrNearest),
								space = as.character(r.o$chr)
								)
		}		
		else if (output == "both" || output=="b")
		{
			debug =0
			if (debug == 0)
			{
			r.o = cbind(as.character(r.o$peak), as.character(r.o$chr),as.numeric(as.character(r.o$peak_start)), 
						as.numeric(as.character(r.o$peak_end)), as.character(r.o$feature),
						as.numeric(as.character(r.o$feature_start)), as.numeric(as.character(r.o$feature_end)),
						as.character(r.o$strand),  as.character(r.o$overlapFeature),
						as.numeric(as.character(r.o$distancetoFeature)),as.character(r.o$fromOverlappingOrNearest),as.numeric(as.character(r.o$shortestDistance))
						)
			colnames(r.o) = c("name","chr", "peakStart", "peakEnd", "feature_id", "start_position", "end_position",
								"strand", "insideFeature", "distancetoFeature", "fromOverlappingOrNearest", "shortestDistance")
			temp = setdiff(paste(r.o[,1], r.o[,5]), paste(r.n[,1], r.n[,5]))
			if (length(temp) >0)
			{
				r.o.only = r.o[paste(r.o[,1], r.o[,5]) %in% temp,]
				r.o.only = matrix(r.o.only, ncol=12)
				colnames(r.o.only) = c("name","chr", "peakStart", "peakEnd", "feature_id", "start_position", "end_position",
								"strand", "insideFeature", "distancetoFeature", "fromOverlappingOrNearest", "shortestDistance")
				r.both = rbind(r.n, r.o.only)
				r.output = RangedData(IRanges(start=as.numeric(as.character(r.both$peakStart)), 
						end=as.numeric(as.character(r.both$peakEnd)),
						names=paste(as.character(r.both$name),as.character(r.both$feature_id))), 
				peak = as.character(r.both$name),
          		strand = as.character(r.both$strand), 
				feature = as.character(r.both$feature_id),
				start_position= as.numeric(as.character(r.both$start_position)),
				end_position=as.numeric(as.character(r.both$end_position)),
				insideFeature=as.character(r.both$insideFeature), 
				distancetoFeature=as.numeric(as.character(r.both$distancetoFeature)),
				shortestDistance = as.numeric(as.character(r.both$shortestDistance)),
				fromOverlappingOrNearest = as.character(r.both$fromOverlappingOrNearest),
				space = as.character(r.both$chr))
			}
			else
			{
				r.output = RangedData(IRanges(start=as.numeric(as.character(r.n$peakStart)), 
						end=as.numeric(as.character(r.n$peakEnd)),
						names=paste(as.character(r.n$name),as.character(r.n$feature_id))),
          		peak = as.character(r.n$name),
				strand = as.character(r.n$strand), 
				feature = as.character(r.n$feature_id),
				start_position= as.numeric(as.character(r.n$start_position)),
				end_position=as.numeric(as.character(r.n$end_position)),
				insideFeature=as.character(r.n$insideFeature), 
				distancetoFeature=as.numeric(as.character(r.n$distancetoFeature)),
				shortestDistance = as.numeric(as.character(r.n$shortestDistance)),
				fromOverlappingOrNearest = as.character(r.n$fromOverlappingOrNearest),
				space = as.character(r.n$chr)
				)
			}
			}
		}
			r.output
		}