BED2RangedData <- function(data.BED, header=FALSE, ...)
{
    if (missing(data.BED)){
        stop("data.BED is required!")
    }
    if (inherits(data.BED, "character")){
        data.BED <- read.delim(data.BED, header=header, ...)
        header <- FALSE
    }
    if ((!is(data.BED, "data.frame")) || dim(data.BED)[2] <3)
    {
        stop("No valid data.BED passed in, which is a data frame or file as BED 
             format with at least 3 fields in the order of: chromosome, start 
             and end. Optional fields are name, score and strand etc. Please 
             refer to http://genome.ucsc.edu/FAQ/FAQformat#format1 for details.  
             If it is a file, please make sure the file path is correct!")
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
        if (any(duplicated(names)))
        {
            names = formatC(1:dim(myPeak)[1], 
                            width=nchar(dim(myPeak)[1]), 
                            flag='0')
        }
    }
    else
    {
        names = formatC(1:dim(myPeak)[1], 
                        width=nchar(dim(myPeak)[1]), 
                        flag='0')
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
        strand = as.character(myPeak[,6])
        strand.allowed.characters=c("1","-1","+","-")
        strand.levels=levels(as.factor(strand))
        if(any(!strand.levels %in% strand.allowed.characters))
        {
            stop("only 1, -1, + and - are acceptable for strand info.
                 Please double check you BED data")
        }
        strand[strand== "1"] = "+"
        strand[strand== "-1"] = "-"
    }
    else
    {
        strand = rep("+", dim(myPeak)[1])
    }
    if(exists("RangedData")){
        RangedData(IRanges(start=as.numeric(as.character(myPeak[,2])), 
                           end=as.numeric(as.character(myPeak[,3])), 
                           names = names), 
                   space = sub("chr", "", as.character(myPeak[,1])), 
                   strand = strand, 
                   score=score)
    }else{
        message("RangedData is dropped, please try ?toGRanges.")
    }
}
