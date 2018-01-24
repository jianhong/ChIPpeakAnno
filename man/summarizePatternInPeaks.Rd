\name{summarizePatternInPeaks}
\alias{summarizePatternInPeaks}
\title{
Output a summary of the occurrence of each pattern in the sequences.
}
\description{
Output a summary of the occurrence of each pattern in the sequences.
}
\usage{
summarizePatternInPeaks(patternFilePath, format = "fasta",skip=0L, 
                        BSgenomeName, peaks, outfile, append = FALSE)
}
\arguments{
  \item{patternFilePath}{
A character vector containing the path to the file to read the patterns from.
}
  \item{format}{
Either "fasta" (the default) or "fastq"
}
 \item{skip}{ Single non-negative integer. The number of records of the pattern
 file to skip before beginning to read in records.
}
  
\item{BSgenomeName}{ BSgenome object. Please refer to available.genomes in
BSgenome package for details
}
  \item{peaks}{
\link[GenomicRanges:GRanges-class]{GRanges} or \link[IRanges:RangedData-class]{RangedData} 
containing the peaks
}
  \item{outfile}{
A character vector containing the path to the file to write the summary output.
}
  \item{append}{
TRUE or FALSE, default FALSE
}
}
\value{
 A data frame with 3 columns as n.peaksWithPattern (number of peaks with the 
 pattern),  n.totalPeaks (total number of peaks in the input) and Pattern 
 (the corresponding pattern). The summary will consider both strand 
 (including reverse complement).
}
\author{
Lihua Julie Zhu
}
\examples{
    peaks = RangedData(IRanges(start=c(100, 500), end=c(300, 600), 
                               names=c("peak1", "peak2")), 
                       space=c("NC_008253", "NC_010468"))
    filepath =system.file("extdata", "examplePattern.fa", 
                          package="ChIPpeakAnno")
    library(BSgenome.Ecoli.NCBI.20080805)
    summarizePatternInPeaks(patternFilePath=filepath, format="fasta", 
                            skip=0L, BSgenomeName=Ecoli, peaks=peaks)
}
\keyword{ misc }