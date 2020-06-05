#' Convert other common IDs to entrez gene ID.
#' 
#' Convert other common IDs such as ensemble gene id, gene symbol, refseq id to
#' entrez gene ID leveraging organism annotation dataset.  For example,
#' org.Hs.eg.db is the dataset from orgs.Hs.eg.db package for human, while
#' org.Mm.eg.db is the dataset from the org.Mm.eg.db package for mouse.
#' 
#' 
#' @param IDs a vector of IDs such as ensembl gene ids
#' @param orgAnn organism annotation dataset such as org.Hs.eg.db
#' @param ID_type type of ID: can be ensemble_gene_id, gene_symbol or refseq_id
#' @return vector of entrez ids
#' @author Lihua Julie Zhu
#' @keywords misc
#' @export
#' @examples
#' 
#' ensemblIDs = c("ENSG00000115956", "ENSG00000071082", "ENSG00000071054",
#'  "ENSG00000115594", "ENSG00000115594", "ENSG00000115598", "ENSG00000170417")
#' library(org.Hs.eg.db)
#' entrezIDs = convert2EntrezID(IDs=ensemblIDs, orgAnn="org.Hs.eg.db",
#'  ID_type="ensembl_gene_id")
#' 
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
        stop("Currently only the following type of IDs are supported:",
             "ensembl_gene_id, refseq_id and gene_symbol!")
    }
    
    if(!is(orgAnn, "AnnDbBimap"))
    {
        stop("orgAnn is not a valid annotation dataset! ",
             "For example, orgs.Hs.eg.db package for human",
			 "and the org.Mm.eg.db package for mouse.")
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
