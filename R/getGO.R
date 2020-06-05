#' Obtain gene ontology (GO) terms for given genes
#' 
#' Obtain gene ontology (GO) terms useing GO gene mapping package such as
#' org.Hs.db.eg to obtain the GO annotation.
#' 
#' 
#' @param all.genes A character vector of feature IDs
#' @param orgAnn Organism annotation package such as org.Hs.eg.db for human and
#' org.Mm.eg.db for mouse, org.Dm.eg.db for fly, org.Rn.eg.db for rat,
#' org.Sc.eg.db for yeast and org.Dr.eg.db for zebrafish
#' @param ID_type The feature type in annotatedPeak such as ensembl_gene_id,
#' refseq_id, gene_symbol
#' @param writeTo File path for output table
#' @return An invisible table with genes and GO terms.
#' @author Lihua Julie Zhu
#' @seealso getEnrichedGO
#' @keywords misc
#' @export
#' @importFrom AnnotationDbi mappedkeys Definition Ontology Term
#' @importFrom utils write.table
#' @examples
#' 
#'   if (interactive()) {
#'      data(annotatedPeak)
#'      library(org.Hs.eg.db)
#'      getGO(annotatedPeak[1:6]$feature, 
#'           orgAnn="org.Hs.eg.db", 
#'           ID_type="ensembl_gene_id")
#' }
#' 
getGO <- function(all.genes, orgAnn = "org.Hs.eg.db", 
                  writeTo, ID_type="gene_symbol"){

  GOgenome = sub(".db", "", orgAnn)
  
  cov2EntrezID <- function(IDs, orgAnn, ID_type = "ensembl_gene_id") {
    
    GOgenome = sub(".db", "", orgAnn)
    
    orgAnn <- switch(ID_type, ensembl_gene_id = "ENSEMBL2EG", 
                     
                     gene_symbol = "SYMBOL2EG", refseq_id = "REFSEQ2EG", 
                     
                     "BAD_NAME")
    
    if (orgAnn == "BAD_NAME") 
      
      stop("Currently only the following type of IDs are supported:\n 
           ensembl_gene_id, refseq_id  and gene_symbol!")
    
    orgAnn <- get(paste(GOgenome, orgAnn, sep = ""))
    
    if (!is(orgAnn, "AnnDbBimap")) {
      
      stop("orgAnn is not a valid annotation dataset! \n 
           For example, orgs.Hs.eg.db package for human and \n 
           the org.Mm.eg.db package for mouse.")
      
    }
    
    IDs <- unique(IDs[!is.na(IDs)])
    
    ids <- mget(IDs, orgAnn, ifnotfound = NA)
    
    ids <- lapply(ids, `[`, 1)
    
    ids <- unlist(ids)
    
    ids <- unique(ids[!is.na(ids)])
    
    
    
    ids
    
  }
  
  
  
  entrezIDs <- cov2EntrezID(as.character(all.genes), orgAnn, ID_type)
  
  
  
  if (length(entrezIDs) < 2) {
    
    stop("The number of gene is less than 2. 
         Please double check your feature_id_type.")
    
  }
  
  goAnn <- get(paste(GOgenome, "GO", sep = ""))
  
  mapped_genes <- mappedkeys(goAnn)
  
  totalN.genes = length(unique(mapped_genes))
  
  thisN.genes = length(unique(entrezIDs))
  
  xx <- mget(mapped_genes, goAnn, ifnotfound = NA)
  
  all.GO <- cbind(matrix(unlist(unlist(xx)), ncol = 3, byrow = TRUE), 
                  
                  rep(names(xx), elementNROWS(xx)))
  
  all.GO <- unique(all.GO)
  
  this.GO <- all.GO[all.GO[, 4] %in% entrezIDs, , drop = FALSE]
  
  
  
  annoTerms <- function(goids) {
    
    if (length(goids) < 1) {
      
      goterm <- matrix(ncol = 4)
      
    }
    
    else {
      
      goids <- as.character(goids)
      
      terms <- Term(goids)
      
      definition <- Definition(goids)
      
      ontology <- Ontology(goids)
      
      goterm <- cbind(goids, terms[match(goids, names(terms))], 
                      
                      definition[match(goids, names(definition))]) 
      
      #, ontology[match(goids, names(ontology))])
      
    }
    
    colnames(goterm) <- c("go.id", "go.term", "Definition" 
                          
    )
    
    rownames(goterm) <- NULL
    
    goterm
    
  }
  
  
  
  this.GO.withTerms <- unique(annoTerms(this.GO[,1]))
  
  colnames(this.GO) = c("go.id", "Evidence", "Ontology", "EntrezID")
  
  
  
  IDs <- as.character(all.genes)
  
  GOgenome = sub(".db", "", orgAnn)
  
  orgAnn2 <- switch(ID_type, ensembl_gene_id = "ENSEMBL2EG", 
                    
                    gene_symbol = "SYMBOL2EG", refseq_id = "REFSEQ2EG", 
                    
                    "BAD_NAME")
  
  if (orgAnn2 == "BAD_NAME") 
    
    stop("Currently only the following type of IDs are supported:\n 
         ensembl_gene_id, refseq_id and gene_symbol!")
  
  orgAnn2 <- get(paste(GOgenome, orgAnn2, sep = ""))
  
  if (!is(orgAnn2, "AnnDbBimap")) {
    
    stop("orgAnn is not a valid annotation dataset! \n 
         For example, orgs.Hs.eg.db package for human and \n 
         the org.Mm.eg.db package for mouse.")
    
  }
  
  
  
  ids <- mget(IDs, orgAnn2, ifnotfound = NA)
  
  ids.withInputIDs <- cbind(names(unlist(ids)), unlist(ids))
  
  ids.withInputIDs <- subset(ids.withInputIDs, !is.na(ids.withInputIDs[,2]))
  
  
  
  this.GO.withTermsEntrezIDs <- merge(this.GO, this.GO.withTerms)
  
  #dim(this.GO.withTermsEntrezIDs)
  
  
  
  this.GO.withTermsEntrezIDs <- 
    cbind(ids.withInputIDs[match(this.GO.withTermsEntrezIDs[, 4], 
                                 ids.withInputIDs[,2]), 1],
          this.GO.withTermsEntrezIDs)
  
  colnames(this.GO.withTermsEntrezIDs)[1] = "Alternate.ID"
  
  
  
  
  
  if(!missing(writeTo)) write.table(this.GO.withTermsEntrezIDs,
                                    file =writeTo, sep="\t",
                                    row.names = FALSE)
  
  return(invisible(this.GO.withTermsEntrezIDs))
}
