\name{summarizePatternInPeaks}
\alias{summarizePatternInPeaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Output a summary of the occurrence of each pattern in the sequences.
}
\description{
Output a summary of the occurrence of each pattern in the sequences.
}
\usage{
summarizePatternInPeaks(patternFilePath, format = "fasta",skip=0L, BSgenomeName, peaks, outfile, append = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{patternFilePath}{
A character vector containing the path to the file to read the patterns from.
}
  \item{format}{
Either "fasta" (the default) or "fastq"
}
 \item{skip}{ Single non-negative integer. The number of records of the pattern file to skip before beginning to read in records.
}
  
\item{BSgenomeName}{ BSgenome object. Please refer to available.genomes in BSgenome package for details
}
  \item{peaks}{
RangedData containing the peaks
}
  \item{outfile}{
A character vector containing the path to the file to write the summary output.
}
  \item{append}{
TRUE or FALSE, default FALSE
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 A data frame with 3 columns as n.peaksWithPattern (number of peaks with the pattern),  n.totalPeaks (total number of peaks in the input) and Pattern (the corresponding pattern).
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Lihua Julie Zhu
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
	peaks = RangedData(IRanges(start=c(100, 500), end=c(300, 600), names=c("peak1", "peak2")), space=c("NC_008253", "NC_010468"))
	filepath =system.file("extdata", "examplePattern.fa", package="ChIPpeakAnno")
	library(BSgenome.Ecoli.NCBI.20080805)
	summarizePatternInPeaks(patternFilePath=filepath, format="fasta", skip=0L, BSgenomeName=Ecoli, peaks=peaks)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }