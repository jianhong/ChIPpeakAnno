GFF2RangedData <- function(data.GFF,header=FALSE)
{
	if (missing(data.GFF) || (class(data.GFF) != "data.frame") || dim(data.GFF)[2] != 9)
	{
		stop("No valid data.GFF passed in, which is a data frame as GFF format with 9 required fields. Please refer to http://genome.ucsc.edu/FAQ/FAQformat#format3 for details.")
	}
	if (header == TRUE)
	{	
		myPeak = data.GFF[-1,]
	}
	else
	{
		myPeak = data.GFF
	}	
	strand = as.character(myPeak[,7])
	strand[strand== "+"] = 1
	strand[strand=="-"] = -1
	names = formatC(1:dim(myPeak)[1], width=nchar(dim(myPeak)[1]), flag='0')
	
	RangedData(IRanges(start=as.numeric(as.character(myPeak[,4])), end=as.numeric(as.character(myPeak[,5])), names = names), space = sub("chr", "", as.character(myPeak[,1])), strand = strand, score=as.numeric(as.character(myPeak[,6])))
}
