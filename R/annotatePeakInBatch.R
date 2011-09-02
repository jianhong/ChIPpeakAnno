annotatePeakInBatch <-
function (myPeakList, mart, featureType = c("TSS", "miRNA", "Exon"), 
    AnnotationData, output = c("nearestStart", "overlapping", 
        "both"), multiple = c(TRUE,FALSE), maxgap = 0, PeakLocForDistance=c("start","middle","end"), FeatureLocForDistance=c("TSS","middle","start","end", "geneEnd"), select=c("all", "first", "last", "arbitrary")) 
{
    featureType = match.arg(featureType)
    if (missing(PeakLocForDistance))
    {
	PeakLocForDistance = "start"
    }
  if (missing(FeatureLocForDistance))
    {
	FeatureLocForDistance = "TSS"
    }

    if (missing(output)) {
        output = "nearestStart"
    }
	select= match.arg(select)
    if ((output == "overlapping" || output == "both") && select =="all" && multiple==FALSE) {
		warning("Please use select instead of multiple!") 
		select = "first"
	}
    if (missing(myPeakList)) {
        stop("Missing required argument myPeakList!")
    }
    if (class(myPeakList) != "RangedData") {
        stop("No valid myPeakList passed in. It needs to be RangedData object")
    }
    if (missing(AnnotationData)) {
        message("No AnnotationData as RangedData is passed in, so now querying biomart database for AnnotationData ....")
        if (missing(mart) || class(mart) != "Mart") {
            stop("Error in querying biomart database. No valid mart object is passed in! Suggest call getAnnotation before calling annotatePeakInBatch")
        }
        AnnotationData <- getAnnotation(mart, feature = featureType)
        message("Done querying biomart database, start annotating ....Better way would be calling getAnnotation before calling annotatePeakInBatch")
    }
    if (class(AnnotationData) != "RangedData") {
        stop("AnnotationData needs to be RangedData object")
    }
    TSS.ordered <- AnnotationData
    rm(AnnotationData)
    if (!length(rownames(TSS.ordered))) {
        rownames(TSS.ordered) = formatC(1:dim(TSS.ordered)[1], 
            width = nchar(dim(TSS.ordered)[1]), flag = "0")
    }
if (length(TSS.ordered$strand) == length(start(TSS.ordered))) 
{ 
	r2 = cbind(rownames(TSS.ordered), start(TSS.ordered), end(TSS.ordered), as.character(TSS.ordered$strand))
}
else
{
	TSS.ordered$strand = rep("+",length(start(TSS.ordered)))
	r2 = cbind(rownames(TSS.ordered), start(TSS.ordered), end(TSS.ordered), rep("+",length(start(TSS.ordered))))
}
colnames(r2) = c("feature_id", "start_position", "end_position", 
        "strand")
if (FeatureLocForDistance == "middle" || FeatureLocForDistance == "m")
{
	FeatureLoc = unlist(lapply(1:dim(r2)[1], function(i)
		{
		round(mean(c(as.numeric(r2[i,2]),as.numeric(r2[i,3]))))
	}))
	FeatureLocForDistance = "middle" 
}
else if (FeatureLocForDistance == "start" || FeatureLocForDistance == "s")
{
	FeatureLoc = r2[,2]
	FeatureLocForDistance = "start"
}
else if (FeatureLocForDistance == "end" || FeatureLocForDistance == "e")
{
	FeatureLoc =r2[,3]
	FeatureLocForDistance = "end"
}
else if (FeatureLocForDistance == "geneEnd" || FeatureLocForDistance == "g")
{
	FeatureLoc = unlist(lapply(1:dim(r2)[1], function(i) 
		{
			 if (as.character(r2[i,4]) == '+' || as.character(r2[i,4])== '1' || as.character(r2[i,4]) == '*' ) {
				r2[i,3]
			} else {
				r2[i,2]
			}
		}))
	FeatureLocForDistance = "geneEnd"
}
else
{
	FeatureLoc = unlist(lapply(1:dim(r2)[1], function(i) 
		{
			 if (as.character(r2[i,4]) == '+' || as.character(r2[i,4])== '1' || as.character(r2[i,4]) == '*' ) {
				r2[i,2]
			} else {
				r2[i,3]
			}
		}))
		FeatureLocForDistance = "TSS"
 }
r2 = cbind(r2, FeatureLoc)

allChr.Anno = unique(as.character(IRanges:::space(TSS.ordered)))
    numberOfChromosome = length(unique(IRanges:::space(myPeakList)))
    if (!length(rownames(myPeakList))) {
        rownames(myPeakList) = formatC(1:dim(myPeakList)[1], 
            width = nchar(dim(myPeakList)[1]), flag = "0")
    }
    allChr = unique(as.character(IRanges:::space(myPeakList)))
    allChr = sub(" +", "", allChr)
    if (length(grep("chr", allChr, fixed = TRUE)) > 0 && length(grep("chr", 
        allChr.Anno, fixed = TRUE)) == 0) {
        allChr = sub("chr", "", allChr)
        myPeakList = RangedData(IRanges(start = start(myPeakList), 
            end = end(myPeakList), names = rownames(myPeakList)), 
            space = sub("chr", "", IRanges:::space(myPeakList)))
    }
    if (length(grep("chr", allChr, fixed = TRUE)) == 0 && length(grep("chr", 
        allChr.Anno, fixed = TRUE)) > 0) {
        allChr = paste("chr", allChr, sep = "")
        myPeakList = RangedData(IRanges(start = start(myPeakList), 
            end = end(myPeakList), names = rownames(myPeakList)), 
            space = paste("chr",IRanges:::space(myPeakList), sep = ""))
    }
    if (output == "nearestStart" || output == "both" || output =="n" || output == "b") {
        z1 = cbind(as.character(rownames(myPeakList)), as.character(IRanges:::space(myPeakList)), 
            start(myPeakList), end(myPeakList))
        colnames(z1) = c("name", "chr", "peakStart", "peakEnd")
        z1[, 2] = sub(" +", "", z1[, 2])

if (PeakLocForDistance == "middle" || PeakLocForDistance == "m")
{
	PeakLoc = unlist(lapply(1:dim(z1)[1], function(i)
		{
		round(mean(c(as.numeric(as.character(z1[i,3])),as.numeric(as.character(z1[i,4])))))
	}))
}
else if (PeakLocForDistance == "start" || PeakLocForDistance == "s")
{
	PeakLoc = as.numeric(as.character(z1[,3]))
}
else if (PeakLocForDistance == "end" || PeakLocForDistance == "e")
{
	PeakLoc =as.numeric(as.character(z1[,4]))
}
else
{
	PeakLoc = as.numeric(as.character(z1[,3]))
}
z1 = cbind(z1, PeakLoc)
TSS.ordered$FeatureLoc = FeatureLoc
myPeakList$PeakLoc = PeakLoc

        plusAnno = TSS.ordered[as.character(TSS.ordered$strand) %in% c("+", "*", "1") , ]
        minusAnno = TSS.ordered[as.character(TSS.ordered$strand) %in% c("-1","-"), ]
        r1 = do.call(rbind, lapply(seq_len(numberOfChromosome), 
            function(i) {
                chr = allChr[i]
                if (chr %in% allChr.Anno) {
                  featureStart = as.numeric(TSS.ordered[chr]$FeatureLoc)
                  peakLoc = as.numeric(myPeakList[chr]$PeakLoc)
	       peakStart = as.numeric(start(myPeakList[chr]))
                  peakEnd = as.numeric(end(myPeakList[chr]))
                  name = rownames(myPeakList[chr])
                  peakRanges = IRanges(start = peakStart, end = peakEnd, 
                    names = name)
                  featureID = rownames(TSS.ordered[chr])
                  featureRanges = IRanges(start = featureStart, 
                    end = featureStart, names = featureID)
                  nearestFeature = featureRanges[nearest(peakRanges, 
                    featureRanges)]
# data.frame(name = name, chr = rep(chr, length(peakStart)), 
#     PeakLoc = peakLoc, peakStart = peakStart, peakEnd = peakEnd, 
                   data.frame(name=name, feature_id = names(nearestFeature))
                }
            }))
        if (length(r1) > 0) {
            r3 = merge(r1, r2, by = "feature_id")
            r = merge(r3, z1, all.y = TRUE)
            r.n = r
	if (FeatureLocForDistance =="TSS" || FeatureLocForDistance =="geneEnd")
	{
            	r11 = r[!is.na(r$strand) & (as.character(r$strand) == "1" | as.character(r$strand) == 
                "+" | as.character(r$strand) == 
                "*"), ]
            	r22 = r[!is.na(r$strand) & (as.character(r$strand) == "-1" | as.character(r$strand) == 
                "-"), ]
          	}
	else
	{
		r11 = r[!is.na(r$strand),]
	}
	r33 = r[is.na(r$strand), ]
           r33$insideFeature = replicate(length(r33$name), NA)
           r33$distancetoFeature = replicate(length(r33$name), NA)

            if (dim(r11)[1] > 0) {
                distancetoFeature = as.numeric(as.character(r11$PeakLoc)) - 
                  as.numeric(as.character(r11$FeatureLoc))
                length = as.numeric(as.character(r11$end_position)) - 
                  as.numeric(as.character(r11$start_position))
                insideFeature = do.call(c, lapply(seq(from = 1, 
                  to = dim(r11)[1], by = 1), function(i) {
                  if (as.numeric(as.character(r11$peakStart[i])) >= 
                    as.numeric(as.character(r11$start_position[i])) && 
                    as.numeric(as.character(r11$peakStart[i])) <= 
                      as.numeric(as.character(r11$end_position[i]))) {
                    if (as.numeric(as.character(r11$peakEnd[i])) >= 
                      as.numeric(as.character(r11$start_position[i])) && 
                      as.numeric(as.character(r11$peakEnd[i])) <= 
                        as.numeric(as.character(r11$end_position[i]))) {
                      "inside"
                    }
                    else {
                      "overlapEnd"
                    }
                  }
                  else if (as.numeric(as.character(r11$peakEnd[i])) >= 
                    as.numeric(as.character(r11$start_position[i])) && 
                    as.numeric(as.character(r11$peakEnd[i])) <= 
                      as.numeric(as.character(r11$end_position[i]))) {
                    "overlapStart"
                  }
                  else if (as.numeric(as.character(r11$peakEnd[i])) >= 
                    as.numeric(as.character(r11$end_position[i])) && 
                    as.numeric(as.character(r11$peakStart[i])) <= 
                      as.numeric(as.character(r11$start_position[i]))) {
                    "includeFeature"
                  }
                  else if (as.numeric(as.character(r11$peakEnd[i])) < 
                    as.numeric(as.character(r11$start_position[i]))) {
                    "upstream"
                  }
                  else if (as.numeric(as.character(r11$peakStart[i])) > 
                    as.numeric(as.character(r11$end_position[i]))) {
                    "downstream"
                  }
                }))
                r11$insideFeature = insideFeature
                r11$distancetoFeature = distancetoFeature
                r.n = r11
            }
            if ((FeatureLocForDistance =="TSS" || FeatureLocForDistance =="geneEnd") && dim(r22)[1] > 0 ){
		distancetoFeature = as.numeric(as.character(r22$FeatureLoc)) - 
                  as.numeric(as.character(r22$PeakLoc))
	length = as.numeric(as.character(r22$end_position)) - 
                  as.numeric(as.character(r22$start_position))
                insideFeature = do.call(c, lapply(seq(from = 1, 
                  to = dim(r22)[1], by = 1), function(i) {
                  if (as.numeric(as.character(r22$peakStart[i])) >= 
                    as.numeric(as.character(r22$start_position[i])) && 
                    as.numeric(as.character(r22$peakStart[i])) <= 
                      as.numeric(as.character(r22$end_position[i]))) {
                    if (as.numeric(as.character(r22$peakEnd[i])) >= 
                      as.numeric(as.character(r22$start_position[i])) && 
                      as.numeric(as.character(r22$peakEnd[i])) <= 
                        as.numeric(as.character(r22$end_position[i]))) {
                      "inside"
                    }
                    else {
                      "overlapStart"
                    }
                  }
                  else if (as.numeric(as.character(r22$peakEnd[i])) >= 
                    as.numeric(as.character(r22$start_position[i])) && 
                    as.numeric(as.character(r22$peakEnd[i])) <= 
                      as.numeric(as.character(r22$end_position[i]))) {
                    "overlapEnd"
                  }
                  else if (as.numeric(as.character(r22$peakEnd[i])) >= 
                    as.numeric(as.character(r22$end_position[i])) && 
                    as.numeric(as.character(r22$peakStart[i])) <= 
                      as.numeric(as.character(r22$start_position[i]))) {
                    "includeFeature"
                  }
                  else if (as.numeric(as.character(r22$peakEnd[i])) < 
                    as.numeric(as.character(r22$start_position[i]))) {
                    "downstream"
                  }
                  else if (as.numeric(as.character(r22$peakStart[i])) > 
                    as.numeric(as.character(r22$end_position[i]))) {
                    "upstream"
                  }
                }))
                r22$insideFeature = insideFeature
                r22$distancetoFeature = distancetoFeature
                if (dim(r11)[1] > 0) {
                  r.n = rbind(r.n, r22)
                }
                else {
                  r.n = r22
                }
            }
            if (dim(r33)[1] > 0) {
                if (dim(r11)[1] > 0 || dim(r22)[1] > 0) {
                  r.n = rbind(r.n, r33)
                }
                else {
                  r.n = r33
                }
            }
            r.n$fromOverlappingOrNearest = rep("NearestStart", 
                dim(r.n)[1])
            shortestDistance = apply(cbind(abs(as.numeric(as.character(r.n$start_position)) - 
                as.numeric(as.character(r.n$peakEnd))), abs(as.numeric(as.character(r.n$start_position)) - 
                as.numeric(as.character(r.n$peakStart))), abs(as.numeric(as.character(r.n$end_position)) - 
                as.numeric(as.character(r.n$peakEnd))), abs(as.numeric(as.character(r.n$end_position)) - 
                as.numeric(as.character(r.n$peakStart)))), 1, 
                min)
            r.n$shortestDistance = shortestDistance
	   strand = do.call(c, lapply(seq_len(dim(r.n)[1]), 
                function(i) {
                  if (is.na(r.n$strand[i]) || as.character(r.n$strand[i]) == "1" || as.character(r.n$strand[i]) == 
                    "+" ||  as.character(r.n$strand[i]) == 
                "*") {
                    "+"
                  }
                  else {
                   "-"
                  }
                }))
            r.n$strand = strand
        }
    }
    if (output == "overlapping" || output == "both" || output =="o" || output == "b") {
        r.o = findOverlappingPeaks(myPeakList, TSS.ordered, maxgap = maxgap, 
            select=select, NameOfPeaks1 = "peak", NameOfPeaks2 = "feature")$OverlappingPeaks
        if (dim(r.o)[1] > 0) {
            r.o$fromOverlappingOrNearest = rep("Overlapping", 
                dim(r.o)[1])
if (PeakLocForDistance == "middle" || PeakLocForDistance == "m")
{
	PeakLoc = unlist(lapply(1:dim(r.o)[1], function(i)
		{
		round(mean(c(as.numeric(as.character(r.o[i,7])),as.numeric(as.character(r.o[i,8])))))
	}))
}
else if (PeakLocForDistance == "start" || PeakLocForDistance == "s")
{
	PeakLoc = as.numeric(as.character(r.o[,7]))
}
else if (PeakLocForDistance == "end" || PeakLocForDistance == "e")
{
	PeakLoc = as.numeric(as.character(r.o [,8]))
}
else
{
	PeakLoc = as.numeric(as.character(r.o[,7]))
}

if (FeatureLocForDistance == "middle" || FeatureLocForDistance == "m")
{
	FeatureLoc = unlist(lapply(1:dim(r.o)[1], function(i)
		{
		round(mean(c(as.numeric(as.character(r.o[i,4])),as.numeric(as.character(r.o[i,5])))))
	}))
}
else if (FeatureLocForDistance == "start" || FeatureLocForDistance == "s")
{
	FeatureLoc =  as.numeric(as.character(r.o[,4]))
}
else if (FeatureLocForDistance == "end" || FeatureLocForDistance == "e")
{
	FeatureLoc = as.numeric(as.character(r.o[,5]))
}
else if (FeatureLocForDistance == "geneEnd" || FeatureLocForDistance == "g")
{
	FeatureLoc = unlist(lapply(1:dim(r2)[1], function(i) 
		{
			 if (as.character(r2[i,4]) == '+' || as.character(r2[i,4])== '1' || as.character(r2[i,4]) == '*' ) {
				r2[i,3]
			} else {
				r2[i,2]
			}
		}))
}
else
{
	FeatureLoc = unlist(lapply(1:dim(r.o)[1], function(i) 
		{
			 if (as.character(r.o[i,6]) == '+' || as.character(r.o[i,6])== '1' || as.character(r.o[i,6]) == '*') {
				 as.numeric(as.character(r.o[i,4]))
			} else {
				 as.numeric(as.character(r.o[i,5]))
			}
		}))
 }
r.o$FeatureLoc = FeatureLoc
r.o$PeakLoc = PeakLoc
#TSS.ordered$FeatureLoc = FeatureLoc
#myPeakList$PeakLoc = PeakLoc

 strand = do.call(c, lapply(seq_len(dim(r.o)[1]), 
                function(i) {
                  if (is.na(r.o$strand[i]) || as.character(r.o$strand[i]) == "1" || as.character(r.o$strand[i]) == 
                    "+" ||  as.character(r.o$strand[i]) == 
                "*") {
                    "+"
                  }
                  else {
                   "-"
                  }
                }))
            r.o$strand = strand


            distancetoFeature = do.call(c, lapply(seq_len(dim(r.o)[1]), 
                function(i) {
                  if (as.character(r.o$strand[i]) == "1" || as.character(r.o$strand[i]) == 
                    "+" ||  as.character(r.o$strand[i]) == 
                "*" || (FeatureLocForDistance != "TSS" && FeatureLocForDistance != "geneEnd")) {
                    as.numeric(as.character(r.o$PeakLoc[i])) - 
                      as.numeric(as.character(r.o$FeatureLoc[i]))
                  }
                  else {
                    as.numeric(as.character(r.o$FeatureLoc[i])) - 
                      as.numeric(as.character(r.o$PeakLoc[i]))
                  }
                }))
            r.o$distancetoFeature = distancetoFeature
        }
    }
    if (output == "nearestStart" || output == "n" || ((output == "both" || output == "b") && dim(r.o)[1] == 0)) {
        if (length(r1) > 0) {
		r.n = as.data.frame(r.n)
            r.output = RangedData(IRanges(start = as.numeric(as.character(r.n$peakStart)), 
                end = as.numeric(as.character(r.n$peakEnd)), 
                names = paste(as.character(r.n$name), as.character(r.n$feature_id))), 
                peak = as.character(r.n$name), strand = as.character(r.n$strand), 
                feature = as.character(r.n$feature_id), start_position = as.numeric(as.character(r.n$start_position)), 
                end_position = as.numeric(as.character(r.n$end_position)), 
                insideFeature = as.character(r.n$insideFeature), 
                distancetoFeature = as.numeric(as.character(r.n$distancetoFeature)), 
                shortestDistance = as.numeric(as.character(r.n$shortestDistance)), 
                fromOverlappingOrNearest = as.character(r.n$fromOverlappingOrNearest), 
                space = as.character(r.n$chr))
        }
        else {
            r.output = myPeakList
        }
    }
    else if ((output == "overlapping" || output == "o") || ((output == "both" || output == "b") && length(r1) == 0)) {
        if (dim(r.o)[1] > 0) {
            r.output = RangedData(IRanges(start = as.numeric(as.character(r.o$peak_start)), 
                end = as.numeric(as.character(r.o$peak_end)), 
                names = paste(as.character(r.o$peak), as.character(r.o$feature))), 
                peak = as.character(r.o$peak), strand = as.character(r.o$strand), 
                feature = as.character(r.o$feature), start_position = as.numeric(as.character(r.o$feature_start)), 
                end_position = as.numeric(as.character(r.o$feature_end)), 
                insideFeature = as.character(r.o$overlapFeature), 
                distancetoFeature = as.numeric(as.character(r.o$distancetoFeature)), 
                shortestDistance = as.numeric(as.character(r.o$shortestDistance)), 
                fromOverlappingOrNearest = as.character(r.o$fromOverlappingOrNearest), 
                space = as.character(r.o$chr))
        }
        else {
            r.output = myPeakList
        }
    }
    else if (output == "both" || output == "b") {
	r.n = as.data.frame(r.n)
	r.n1 = r.n
        debug = 0
        if (debug == 0) {
	   
	    r.n = cbind(as.character(r.n$name),as.character(r.n$chr),  
			as.numeric(as.character(r.n$peakStart)), as.numeric(as.character(r.n$peakEnd)), 
                	as.character(r.n$feature_id), as.numeric(as.character(r.n$start_position)), 
                	as.numeric(as.character(r.n$end_position)), as.character(r.n$strand), 
                	as.character(r.n$insideFeature), as.numeric(as.character(r.n$distancetoFeature)), 
                	as.character(r.n$fromOverlappingOrNearest), as.numeric(as.character(r.n$shortestDistance)))
	  colnames(r.n) = c("name", "chr", "peakStart", "peakEnd", 
                "feature_id", "start_position", "end_position", 
                "strand", "insideFeature", "distancetoFeature", 
                "fromOverlappingOrNearest", "shortestDistance")
            r.o = cbind(as.character(r.o$peak), as.character(r.o$chr), 
                as.numeric(as.character(r.o$peak_start)), as.numeric(as.character(r.o$peak_end)), 
                as.character(r.o$feature), as.numeric(as.character(r.o$feature_start)), 
                as.numeric(as.character(r.o$feature_end)), as.character(r.o$strand), 
                as.character(r.o$overlapFeature), as.numeric(as.character(r.o$distancetoFeature)), 
                as.character(r.o$fromOverlappingOrNearest), as.numeric(as.character(r.o$shortestDistance)))
            colnames(r.o) = c("name", "chr", "peakStart", "peakEnd", 
                "feature_id", "start_position", "end_position", 
                "strand", "insideFeature", "distancetoFeature", 
                "fromOverlappingOrNearest", "shortestDistance")
            temp = setdiff(paste(r.o[, 1], r.o[, 5]), paste(r.n[, 
                1], r.n[, 5]))
            if (length(temp) > 0) {
                r.o.only = r.o[paste(r.o[, 1], r.o[, 5]) %in% 
                  temp, ]
                r.o.only = matrix(r.o.only, ncol = 12)
                colnames(r.o.only) = c("name", "chr", "peakStart", 
                  "peakEnd", "feature_id", "start_position", 
                  "end_position", "strand", "insideFeature", 
                  "distancetoFeature", "fromOverlappingOrNearest", 
                  "shortestDistance")
                r.both = as.data.frame(rbind(r.n, r.o.only))
                r.output = RangedData(IRanges(start = as.numeric(as.character(r.both$peakStart)), 
                  end = as.numeric(as.character(r.both$peakEnd)), 
                  names = paste(as.character(r.both$name), as.character(r.both$feature_id))), 
                  peak = as.character(r.both$name), strand = as.character(r.both$strand), 
                  feature = as.character(r.both$feature_id), 
                  start_position = as.numeric(as.character(r.both$start_position)), 
                  end_position = as.numeric(as.character(r.both$end_position)), 
                  insideFeature = as.character(r.both$insideFeature), 
                  distancetoFeature = as.numeric(as.character(r.both$distancetoFeature)), 
                  shortestDistance = as.numeric(as.character(r.both$shortestDistance)), 
                  fromOverlappingOrNearest = as.character(r.both$fromOverlappingOrNearest), 
                  space = as.character(r.both$chr))
            }
            else {
                r.output = RangedData(IRanges(start = as.numeric(as.character(r.n1$peakStart)), 
                  end = as.numeric(as.character(r.n1$peakEnd)), 
                  names = paste(as.character(r.n1$name), as.character(r.n1$feature_id))), 
                  peak = as.character(r.n1$name), strand = as.character(r.n1$strand), 
                  feature = as.character(r.n1$feature_id), start_position = as.numeric(as.character(r.n1$start_position)), 
                  end_position = as.numeric(as.character(r.n1$end_position)), 
                  insideFeature = as.character(r.n1$insideFeature), 
                  distancetoFeature = as.numeric(as.character(r.n1$distancetoFeature)), 
                  shortestDistance = as.numeric(as.character(r.n1$shortestDistance)), 
                  fromOverlappingOrNearest = as.character(r.n1$fromOverlappingOrNearest), 
                  space = as.character(r.n1$chr))
            }
        }
    }
r.output
}
