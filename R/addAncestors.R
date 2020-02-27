addAncestors <- function(go.ids, ontology=c("bp","cc", "mf"))
{
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
        xx <- as.list(GOBPANCESTOR)
    }
    else if (ontology =="cc")
    {
        xx <- as.list(GOCCANCESTOR)
    }
    else if (ontology =="mf")
    {
        xx <- as.list(GOMFANCESTOR)
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