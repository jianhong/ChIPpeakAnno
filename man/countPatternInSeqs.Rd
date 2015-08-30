\name{countPatternInSeqs}
\alias{countPatternInSeqs}
\title{
Output total number of patterns found in the input sequences
}
\description{
Output total number of patterns found in the input sequences
}
\usage{
countPatternInSeqs(pattern, sequences)
}
\arguments{
  \item{pattern}{
DNAstringSet object
}
  \item{sequences}{
 a vector of sequences
}
}
\value{
Total number of occurrence of the pattern in the sequences
}
\author{
Lihua Julie Zhu
}
\seealso{
summarizePatternInPeaks, translatePattern
}
\examples{
  filepath = 
    system.file("extdata", "examplePattern.fa", package="ChIPpeakAnno")
  dict = readDNAStringSet(filepath = filepath, format="fasta", use.names=TRUE)
  sequences = c("ACTGGGGGGGGCCTGGGCCCCCAAAT", 
                "AAAAAACCCCTTTTGGCCATCCCGGGACGGGCCCAT", 
                "ATCGAAAATTTCC")
  countPatternInSeqs(pattern=dict[1], sequences=sequences)
  countPatternInSeqs(pattern=dict[2], sequences=sequences)
  pattern = DNAStringSet("ATNGMAA")
  countPatternInSeqs(pattern=pattern, sequences=sequences)
}
\keyword{ misc }
