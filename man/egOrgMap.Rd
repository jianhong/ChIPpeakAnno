\name{egOrgMap}
\alias{egOrgMap}
\title{
  map organism annotation dataset to specie name or revese. 
}
\description{
  Give a specie name and return the organism annotation dataset name or 
  give a organism annotation dataset name then return the specie name.
}
\usage{
egOrgMap(name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{
  organism annotation dataset or the specie name.
}
}
\value{
  a object of character
}
\author{
  Jianhong Ou
}
\examples{
  egOrgMap("org.Hs.eg.db")
  egOrgMap("Mus musculus")
}
\keyword{ misc }
