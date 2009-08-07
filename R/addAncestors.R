addAncestors <- function(go.ids, ontology="bp")
{
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
	go.ids.withAncestor = intersect(as.character(go.ids), names(xx))
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
	children = cbind(go.ids[go.ids %in% go.ids.withAncestor])
	colnames(children) = c("child")
	go.all = merge(children, Ancestors, by="child")
	c(as.character(go.ids),as.character(go.all$ancestor))
  }
  else
  {
	as.character(go.ids)
 }
}