\name{getUniqueGOidCount}
\alias{getUniqueGOidCount}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ get the count for each unique GO ID }
\description{
  get the count for each unique GO ID
}
\usage{
getUniqueGOidCount(goList)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{goList}{a set of GO terms as character vector}
}
\value{
  a list with 2 variables
  \item{GOterm}{a vector of GO terms as character vector}
  \item{GOcount}{counts corresponding to the above GOterm as numeric vector}
}
\author{ Lihua Julie Zhu}
\note{ 
internal function not intended to be called directly by users
}
\seealso{ getEnrichedGO }
\examples{
goList= c("GO:0000075", "GO:0000082","GO:0000082","GO:0000122","GO:0000122",
           "GO:0000075","GO:0000082","GO:0000082","GO:0000122","GO:0000122",
           "GO:0000122","GO:0000122","GO:0000075", "GO:0000082","GO:000012")
 
getUniqueGOidCount(goList)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
