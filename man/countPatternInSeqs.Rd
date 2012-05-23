\name{countPatternInSeqs}
\alias{countPatternInSeqs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Output total number of patterns found in the input sequences
}
\description{
Output total number of patterns found in the input sequences
}
\usage{
countPatternInSeqs(pattern, sequences)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pattern}{
DNAstringSet object
}
  \item{sequences}{
 a vector of sequences
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Total number of occurrence of the pattern in the sequences
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
summarizePatternInPeaks, translatePattern
}
\examples{
	filepath = system.file("extdata", "examplePattern.fa", package="ChIPpeakAnno")
	dict = readDNAStringSet(filepath = filepath, format="fasta", use.names=TRUE)
	sequences = c("ACTGGGGGGGGCCTGGGCCCCCAAAT", "AAAAAACCCCTTTTGGCCATCCCGGGACGGGCCCAT", "ATCGAAAATTTCC")
	countPatternInSeqs(pattern=dict[1], sequences=sequences)
	countPatternInSeqs(pattern=dict[2], sequences=sequences)
	pattern = DNAStringSet("ATNGMAA")
	countPatternInSeqs(pattern=pattern, sequences=sequences)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
