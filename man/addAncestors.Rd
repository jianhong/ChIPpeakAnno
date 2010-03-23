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
addAncestors(go.ids, ontology = c("bp", "cc", "mf"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{go.ids}{
 matrix with 4 columns: first column is GO IDs and 4th column is entrez IDs.
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
go.ids = cbind(c("GO:0008150", "GO:0005576", "GO:0003674"),c("ND", "IDA", "ND"), c("BP", "BP", "BP"), c("1", "1", "1"))
addAncestors(go.ids, ontology="bp")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
