GFF2RangedData <- function(data.GFF, header=FALSE, ...)
{
    if (missing(data.GFF)){
        stop("data.GFF is required!")
    }
    if (inherits(data.GFF, "character")){
        data.GFF <- read.delim(data.GFF, header=header, ...)
        header <- FALSE
    }
    if ((class(data.GFF) != "data.frame") || dim(data.GFF)[2] != 9)
    {
        stop("No valid data.GFF passed in, which is a data frame or file as GFF
             format with 9 required fields. Please refer to 
             http://genome.ucsc.edu/FAQ/FAQformat#format3 for details. 
             If it is a file, please make sure the file path is correct!")
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
    strand.allowed.characters=c("1","-1","+","-")
    strand.levels=levels(as.factor(strand))
    if(any(!strand.levels %in% strand.allowed.characters))
    {
        stop("only 1, -1, +, - and . are acceptable for strand info. 
             Please double check you GFF data. If strand entries contain '.', 
             please convert it into '+' or '-'.")
    }
    strand[strand== "1"] = "+"
    strand[strand=="-1"] = "-"
    names = formatC(1:dim(myPeak)[1], width=nchar(dim(myPeak)[1]), flag='0')
    
    if(exists("RangedData")){
        RangedData(IRanges(start=as.numeric(as.character(myPeak[,4])),
                           end=as.numeric(as.character(myPeak[,5])), 
                           names = names), 
                   space = sub("chr", "", as.character(myPeak[,1])), 
                   strand = strand, 
                   score=as.numeric(as.character(myPeak[,6])))
    }else{
        message("RangedData is dropped, please try ?toGRanges.")
    }
}
