#' Add GO IDs of the ancestors for a given vector of GO ids
#' @description Add GO IDs of the ancestors for a given vector of GO IDs 
#' leveraging GO.db 
#' @param go.ids A matrix with 4 columns: first column is GO IDs and 4th 
#' column is entrez IDs.
#' @param ontology bp for biological process, cc for cellular component and 
#' mf for molecular function.
#' @return A vector of GO IDs containing the input GO IDs with the GO IDs of 
#' their ancestors added.
#' @export
#' @author Lihua Julie Zhu
#' @examples 
#' go.ids = cbind(c("GO:0008150", "GO:0005576", "GO:0003674"),
#'                c("ND", "IDA", "ND"), 
#'                c("BP", "BP", "BP"), 
#'                c("1", "1", "1"))
#' library(GO.db)
#' addAncestors(go.ids, ontology="bp")
#' @keywords misc
#' 

addAncestors <- function(go.ids, ontology=c("bp","cc", "mf"))
{
    stopifnot("The 'GO.db' package is required"=
                  requireNamespace("GO.db", quietly = TRUE))
    ontology = match.arg(ontology)
    if (missing(go.ids))
    {
        stop("missing required parameter go.ids!")
    }
    if (!is(go.ids, "matrix") | dim(go.ids)[2] <4)
    {
        stop("go.ids need to be a matrix with at least 4 columns, 
            the first column is go id and the second column is the entrez ID!")
    }
    if (ontology =="bp")
    {
        xx <- as.list(GO.db::GOBPANCESTOR)
    }
    else if (ontology =="cc")
    {
        xx <- as.list(GO.db::GOCCANCESTOR)
    }
    else if (ontology =="mf")
    {
        xx <- as.list(GO.db::GOMFANCESTOR)
    }
    xx <- xx[!is.na(xx)]
    
    if (length(xx) <=0)
    {
        stop("No GO Ancestors! Should not have happened!")
    }
    go.ids.withAncestor = intersect(as.character(go.ids[,1]), names(xx))
    if (length(go.ids.withAncestor) >0)
    {
        Ancestors =	do.call(rbind, lapply(go.ids.withAncestor,function(x1)
        {
            r= xx[names(xx)==x1]
            if (length(r) > 0 )
            {
                cbind(rep(x1, length(r[[1]])), r[[1]])
            }
        }))
        Ancestors = Ancestors[Ancestors[,2] !="all",]
        colnames(Ancestors) = c("child", "ancestor")
        children = cbind(go.ids[,1],go.ids[,4])
        colnames(children) = c("child", "entrezID")
        go.all = merge(children, Ancestors, by="child", all.x=TRUE)
        temp = cbind(c(as.character(go.all$child),
                       as.character(go.all$ancestor)), 
                     c(as.character(go.all$entrezID), 
                       as.character(go.all$entrezID)))
        temp = unique(temp)
        if (length(temp) <3)
        {
            unique(children)
        }
        else
        {
            temp[!is.na(temp[,1]) & temp[,1] != "",]
        }
    }
    else
    {
        unique(cbind(as.character(go.ids[,1]), go.ids[,4]))
    }
}
