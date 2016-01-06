\name{xget}
\alias{xget}
\title{
Return the value from a Bimap objects
}
\description{
Search by name for an Bimap object.
}
\usage{
xget(x, envir, mode, ifnotfound=NA, inherits, 
    output=c("all", "first"))
}
\arguments{
  \item{x, envir, mode, ifnotfound, inherits}{
    see \link[AnnotationDbi]{mget}
  }
  \item{output}{return the all or first item for each query}
}
\value{
    a character vector
}
\author{
    Jianhong Ou
}
\seealso{
 See Also as \code{\link[base]{mget}}, \code{\link[AnnotationDbi]{mget}}
}
\examples{
    library(org.Hs.eg.db)
    xget(as.character(1:10), org.Hs.egSYMBOL)
}
\keyword{ misc }
