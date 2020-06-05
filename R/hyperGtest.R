#' hypergeometric test
#' 
#' hypergeometric test with lower.tail = FALSE used by getEnrichedGO
#' 
#' see phyper for details
#' 
#' @param alltermcount a list with two variables: GOterm and GOcount which is
#' GO terms and corresponding counts in the whole genome
#' @param thistermcount a list with two variables: GOterm and GOcount which is
#' GO terms and corresponding counts in the peak list
#' @param totaltermInGenome number of total GO terms in the whole genome
#' @param totaltermInPeakList number of total GO terms in the peak list
#' @return a list with 6 variables \item{thisterm}{GO term}
#' \item{thistermcount}{count of this GO term in the peak list}
#' \item{thistermtotal}{count of this GO term in the whole genome}
#' \item{pvalue}{pvalue of the hypergeometric test}
#' \item{totaltermInPeakList}{number of total GO terms in the peak list}
#' \item{totaltermInGenome}{number of total GO terms in the whole genome}
#' @note internal function not intended to be used directly by users
#' @author Lihua Julie ZHu
#' @seealso phyper, getEnrichedGO
#' @references Johnson, N. L., Kotz, S., and Kemp, A. W. (1992) Univariate
#' Discrete Distributions, Second Edition. New York: Wiley
#' @keywords internal
#' @export
#' @importFrom stats phyper
#' @examples
#' 
#' goList= c("GO:0000075", "GO:0000082","GO:0000082","GO:0000122",
#'           "GO:0000122","GO:0000075","GO:0000082","GO:0000082",
#'           "GO:0000122","GO:0000122","GO:0000122","GO:0000122",
#'           "GO:0000075", "GO:0000082","GO:000012")
#'   
#' alltermcount = list(GOterm=c("GO:0000075", "GO:0000082", "GO:000012", 
#'                              "GO:0000122"), 
#'                     GOcount=c(100, 200, 10, 10))
#' thistermcount = getUniqueGOidCount(goList)
#' totaltermInPeakList = 15
#' totaltermInGenome = 1000
#' hyperGtest(alltermcount,thistermcount, totaltermInGenome, totaltermInPeakList)
#' 
hyperGtest <- function(alltermcount,thistermcount, 
                       totaltermInGenome, totaltermInPeakList){
    pvalue = as.numeric()
    thistermtotal = as.numeric()
    for (i in 1:length(thistermcount$GOterm))
    {
        m = 
            as.numeric(alltermcount$GOcount[alltermcount$GOterm==
                                                thistermcount$GOterm[i]])
        q = as.numeric(thistermcount$GOcount[i])
        n = as.numeric(totaltermInGenome)
        k = as.numeric(totaltermInPeakList)
        ####### k - number of total GO terms in the peak list
        ####### n - number of GO terms total
        ####### q - number of this GO terms in the peak list
        ####### m - number of this GO terms total
        pvalue[i] = phyper(q-1, m, n-m, k, lower.tail = FALSE, log.p = FALSE)
        thistermtotal[i] = m
    }
    list(thisterm=thistermcount$GOterm, 
         thistermcount=thistermcount$GOcount,
         thistermtotal=thistermtotal, 
         pvalue=pvalue, 
         totaltermInPeakList=k, 
         totaltermInGenome=n)
}

