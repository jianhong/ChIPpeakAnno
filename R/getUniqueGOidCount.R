#' get the count for each unique GO ID
#' 
#' get the count for each unique GO ID
#' 
#' 
#' @param goList a set of GO terms as character vector
#' @return a list with 2 variables \item{GOterm}{a vector of GO terms as
#' character vector} \item{GOcount}{counts corresponding to the above GOterm as
#' numeric vector}
#' @note internal function not intended to be called directly by users
#' @author Lihua Julie Zhu
#' @seealso getEnrichedGO
#' @keywords internal
#' @export
#' @examples
#' 
#' goList= c("GO:0000075", "GO:0000082","GO:0000082","GO:0000122","GO:0000122",
#'            "GO:0000075","GO:0000082","GO:0000082","GO:0000122","GO:0000122",
#'            "GO:0000122","GO:0000122","GO:0000075", "GO:0000082","GO:000012")
#'  
#' getUniqueGOidCount(goList)
#' 
getUniqueGOidCount <-
    function(goList)
    {
        x = goList[order(goList)]
        duplicated.go<-duplicated(x)
        unique.go <- unique(x)
        go.count<- numeric()
        
        count <- 1
        j <- 1
        
        for (i in 2:length(duplicated.go))
        {
            
            if (duplicated.go[i] == FALSE)
            {
                go.count[j] <- count # previous GO
                j <- j + 1
                count <- 1
            }
            else
            {
                count <- count + 1
            }
        }
        go.count[j] <- count #last GO
        list(GOterm=unique.go, GOcount=go.count)
    }

