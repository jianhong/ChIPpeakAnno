\name{hyperGtest}
\alias{hyperGtest}
\title{ hypergeometric test }
\description{
  hypergeometric test with lower.tail = FALSE used by getEnrichedGO
}
\usage{
hyperGtest(alltermcount, thistermcount, totaltermInGenome, totaltermInPeakList)
}
\arguments{
  \item{alltermcount}{ a list with two variables: GOterm and GOcount which is 
  GO terms and corresponding counts in the whole genome }
  \item{thistermcount}{ a list with two variables: GOterm and GOcount which is 
  GO terms and corresponding counts in the peak list}
  \item{totaltermInGenome}{ number of total GO terms in the whole genome}
  \item{totaltermInPeakList}{ number of total GO terms in the peak list }
}
\details{
 see phyper for details
}
\value{
  a list with 6 variables
  \item{thisterm}{GO term}
  \item{thistermcount}{count of this GO term in the peak list}
  \item{thistermtotal}{count of this GO term in the whole genome}
  \item{pvalue}{pvalue of the hypergeometric test}
  \item{totaltermInPeakList}{number of total GO terms in the peak list}
  \item{totaltermInGenome}{number of total GO terms in the whole genome}
}
\references{
Johnson, N. L., Kotz, S., and Kemp, A. W. (1992) Univariate Discrete 
Distributions, Second Edition. New York: Wiley }
\author{ Lihua Julie ZHu }
\note{ 
	internal function not intended to be used directly by users
}
\seealso{ phyper, getEnrichedGO }
\examples{
goList= c("GO:0000075", "GO:0000082","GO:0000082","GO:0000122",
          "GO:0000122","GO:0000075","GO:0000082","GO:0000082",
          "GO:0000122","GO:0000122","GO:0000122","GO:0000122",
          "GO:0000075", "GO:0000082","GO:000012")
  
alltermcount = list(GOterm=c("GO:0000075", "GO:0000082", "GO:000012", 
                             "GO:0000122"), 
                    GOcount=c(100, 200, 10, 10))
thistermcount = getUniqueGOidCount(goList)
totaltermInPeakList = 15
totaltermInGenome = 1000
hyperGtest(alltermcount,thistermcount, totaltermInGenome, totaltermInPeakList)
}
\keyword{ internal }
