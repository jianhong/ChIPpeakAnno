BED2RangedData <- function(data.BED,header=FALSE)
{
	if (missing(data.BED) || (class(data.BED) != "data.frame") || dim(data.BED)[2] <3)
	{
		stop("No valid data.BED passed in, which is a data frame as BED format file with at least 3 fields in the order of: chromosome, start and end. Optional fields are name, score and strand etc. Please refer to http://genome.ucsc.edu/FAQ/FAQformat#format1 for details.")
	}
	if (header == TRUE)
	{	
		myPeak = data.BED[-1,]
	}
	else
	{
		myPeak = data.BED
	}	
	if (dim(myPeak)[2] >=4)
	{
		names = as.character(myPeak[,4])
	}
	else
	{
		names = paste(myPeak[,1], myPeak[,2], myPeak[,3])
	}
	if (dim(myPeak)[2] >= 5)
	{		
		score = as.numeric(as.character(myPeak[,5]))
	}
	else
	{
		score = rep(1, dim(myPeak)[1])
	}
	if (dim(myPeak)[2] >= 6)
	{		
		strand = myPeak[,6]
		strand[strand== "+"] = 1
		strand[strand=="-"] = -1
	}
	else
	{
		strand = rep(1, dim(myPeak)[1])
	}
	RangedData(IRanges(start=as.numeric(as.character(myPeak[,2])), end=as.numeric(as.character(myPeak[,3])), names = names), space = sub("chr", "", as.character(myPeak[,1])), strand = strand, score=score )
}
