\name{HOT.spots}
\alias{HOT.spots}
\docType{data}
\title{
High Occupancy of Transcription Related Factors regions
}
\description{
High Occupancy of Transcription Related Factors regions of human (hg19)
}
\usage{data("HOT.spots")}
\format{
  An object of GRangesList
}
\details{
How to generated the data:

  temp <- tempfile()
  
  url <- "http://metatracks.encodenets.gersteinlab.org"
  
  download.file(file.path(url, "HOT_All_merged.tar.gz"), temp)
  
  temp2 <- tempfile()
  
  download.file(file.path(url, "HOT_intergenic_All_merged.tar.gz"), temp2)
  
  untar(temp, exdir=dirname(temp))
  
  untar(temp2, exdir=dirname(temp))
  
  f <- dir(dirname(temp), "bed$")
  
  HOT.spots <- sapply(file.path(dirname(temp), f), toGRanges, format="BED")
  
  names(HOT.spots) <- gsub("_merged.bed", "", f)
  
  HOT.spots <- sapply(HOT.spots, unname)
  
  HOT.spots <- GRangesList(HOT.spots)
  
  save(list="HOT.spots", 
  
       file="data/HOT.spots.rda", 
       
       compress="xz", compression_level=9)
}
\source{
http://metatracks.encodenets.gersteinlab.org/
}
\references{
 Yip KY, Cheng C, Bhardwaj N, Brown JB, Leng J, Kundaje A, Rozowsky J, Birney
 E, Bickel P, Snyder M, Gerstein M. Classification of human genomic regions 
 based on experimentally determined binding sites of more than 100 
 transcription-related factors. Genome Biol. 2012 Sep 26;13(9):R48. 
 doi: 10.1186/gb-2012-13-9-r48. PubMed PMID: 22950945; 
 PubMed Central PMCID: PMC3491392.
}
\examples{
 data(HOT.spots)
 elementLengths(HOT.spots)
}
\keyword{datasets}
