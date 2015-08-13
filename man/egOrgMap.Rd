\name{egOrgMap}
\alias{egOrgMap}
\title{
  Convert between the name of the organism annotation package ("OrgDb") and the name of the organism. 
}
\description{
  Give a species name and return the organism annotation package name or 
  give an organism annotation package name then return the species name.
}
\usage{
egOrgMap(name)
}
\arguments{
  \item{name}{
  The name of the organism annotation package or the species.
}
}
\value{
  A object of character
}
\author{
  Jianhong Ou
}
\examples{
  egOrgMap("org.Hs.eg.db")
  egOrgMap("Mus musculus")
}
\keyword{ misc }
