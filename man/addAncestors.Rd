\name{addAncestors}
\Rdversion{1.1}
\alias{addAncestors}
\title{
Add GO IDs of the ancestors for a given vector of GO ids
}
\description{
Add GO IDs of the ancestors for a given vector of GO IDs leveraging GO.db package 
}
\usage{
addAncestors(go.ids, ontology = c("bp", "cc", "mf"))
}
\arguments{
  \item{go.ids}{
 A matrix with 4 columns: first column is GO IDs and 4th column is entrez IDs.
}
  \item{ontology}{
 bp for biological process, cc for cellular component and mf for molecular function
}
}
\value{
 A vector of GO IDs containing the input GO IDs with the GO IDs of their ancestors added
}
\author{
Lihua Julie Zhu
}
\examples{
go.ids = cbind(c("GO:0008150", "GO:0005576", "GO:0003674"),c("ND", "IDA", "ND"), 
c("BP", "BP", "BP"), c("1", "1", "1"))
addAncestors(go.ids, ontology="bp")
}
\keyword{ misc }
