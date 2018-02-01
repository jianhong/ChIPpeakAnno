convert2EntrezID <- function(IDs, orgAnn, ID_type="ensembl_gene_id")
{
    GOgenome = sub(".db","",orgAnn)
    if (ID_type == "ensembl_gene_id")
    {
        orgAnn <- get(paste(GOgenome,"ENSEMBL2EG", sep=""))
    }
    else if (ID_type == "gene_symbol")
    {
        orgAnn <- get(paste(GOgenome,"SYMBOL2EG", sep=""))
    }
    else if (ID_type == "refseq_id")
    {
        orgAnn <- get(paste(GOgenome,"REFSEQ2EG", sep=""))
    }
    else
    {
        stop("Currently only the following type of IDs are supported: ensembl_gene_id, refseq_id and gene_symbol!")
    }
    
    if(!is(orgAnn, "AnnDbBimap"))
    {
        stop("orgAnn is not a valid annotation dataset! For example, orgs.Hs.eg.db package for human
			 and the org.Mm.eg.db package for mouse.")
    }
    xx <- as.list(orgAnn)
    IDs = unique(IDs[!is.na(IDs)])
    x = do.call(rbind, lapply(IDs,function(x1)
    {
        r= xx[names(xx)==x1]
        if (length(r) >0)
        {
            r[[1]][1]
        }
    }))
    unique(x)
}