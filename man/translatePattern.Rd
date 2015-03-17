\name{translatePattern}
\alias{translatePattern}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 translate pattern from IUPAC Extended Genetic Alphabet to regular expression
}
\description{
   translate pattern containing the IUPAC nucleotide ambiguity codes to regular expression. 
For example,Y->[C|T], R-> [A|G], S-> [G|C], W-> [A|T], 	K-> [T|U|G], M-> [A|C],
 B-> [C|G|T], D-> [A|G|T], H-> [A|C|T], V-> [A|C|G] and N-> [A|C|T|G].
}
\usage{
translatePattern(pattern)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pattern}{
a character vector with the IUPAC nucleotide ambiguity codes
}
}
\value{
a character vector with the pattern represented as regular expression
}
\author{
Lihua Julie Zhu
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
countPatternInSeqs, summarizePatternInPeaks
}
\examples{
	pattern1 = "AACCNWMK"
	translatePattern(pattern1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
