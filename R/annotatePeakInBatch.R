annotatePeakInBatch <-
function(myPeakList, mart,featureType=c("TSS","miRNA", "Exon"), AnnotationData)
{
		featureType = match.arg(featureType)
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
		
		r2 = cbind(rownames(TSS.ordered), start(TSS.ordered), end(TSS.ordered), TSS.ordered$strand)
		colnames(r2) = c("feature_id", "start_position", "end_position", "strand")
		allChr.Anno = unique(space(TSS.ordered))
		
		#allChr.Anno = sub("chr", '', allChr.Anno)
		
		numberOfChromosome = length(unique(space(myPeakList)))
	
		if (!length(rownames(myPeakList)))
		{
			rownames(myPeakList) = paste(as.character(space(myPeakList)),start(myPeakList),sep=":")
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
		
		r3 = merge(r1,r2, by="feature_id")
		
		r = merge(r3, z1, all.y=TRUE)
		
		r11 = r[!is.na(r$strand) & (r$strand==1 | r$strand =="+"),]
		r22 = r[!is.na(r$strand) & (r$strand==-1 | r$strand =="-"),]
		r33 = r[is.na(r$strand),]
		r33$insideFeature = replicate(length(r33$name), NA)
		r33$distancetoFeature = replicate(length(r33$name), NA)
		
		distancetoFeature =  as.numeric(as.character(r11$peakStart)) - as.numeric(as.character(r11$start_position))
		length =  as.numeric(as.character(r11$end_position)) -  as.numeric(as.character(r11$start_position))
		insideFeature = (distancetoFeature>=0 & distancetoFeature <= length)		
		r11$insideFeature = insideFeature		
		r11$distancetoFeature = distancetoFeature
		
		distancetoFeature =  as.numeric(as.character(r22$end_position)) - as.numeric(as.character(r22$peakStart))
		length =  as.numeric(as.character(r22$end_position)) -  as.numeric(as.character(r22$start_position))
		insideFeature = (distancetoFeature>=0 & distancetoFeature <= length)		
		r22$insideFeature = insideFeature
		r22$distancetoFeature = distancetoFeature
		
		RangedData(IRanges(start=c(as.numeric(as.character(r11$peakStart)),as.numeric(as.character(r22$peakStart)),as.numeric(as.character(r33$peakStart))), 
						end=c(as.numeric(as.character(r11$peakEnd)),as.numeric(as.character(r22$peakEnd)),as.numeric(as.character(r33$peakEnd))),
						names=c(as.character(r11$name),as.character(r22$name),as.character(r33$name))),
          		strand = c(as.character(r11$strand),as.character(r22$strand), r33$strand), 
				feature = c(as.character(r11$feature_id),as.character(r22$feature_id), r33$feature_id),
				start_position= c(as.numeric(as.character(r11$start_position)),as.numeric(as.character(r22$start_position)),r33$start_position),
				end_position=c(as.numeric(as.character(r11$end_position)),as.numeric(as.character(r22$end_position)), r33$end_position),
				insideFeature=c(unlist(r11$insideFeature),unlist(r22$insideFeature),unlist(r33$insideFeature)), 
				distancetoFeature=c(unlist(as.numeric(as.character(r11$distancetoFeature))), unlist(as.numeric(as.character(r22$distancetoFeature))), unlist(r33$distancetoFeature)),
				space = c(as.character(r11$chr),as.character(r22$chr), as.character(r33$chr)))
		}
		