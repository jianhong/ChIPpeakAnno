\name{addAncestors}
\Rdversion{1.1}
\alias{addAncestors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add GO ids of the ancestors for a given vector of GO ids
}
\description{
Add GO ids of the ancestors for a given vector of GO ids leveraging GO.db package 
}
\usage{
addAncestors(go.ids, ontology = "bp")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{go.ids}{
 a vector of GO IDs
}
  \item{ontology}{
 bp for biological process, cc for cellular component and mf for molecular function
}
}
\details{

}
\value{
 a vector of GO IDs containing the input GO IDs with the GO IDs of their ancestors added
}
\references{
}
\author{
Lihua Julie Zhu
}
\note{

}

\seealso{

}
\examples{
go.ids = c("GO:0000001", "GO:0000002", "GO:0000003")
addAncestors(go.ids, ontology="bp")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
