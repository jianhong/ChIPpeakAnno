#' Obtain enriched gene ontology (GO) terms that near the peaks
#' 
#' Obtain enriched gene ontology (GO) terms based on the features near the
#' enriched peaks using GO.db package and GO gene mapping package such as
#' org.Hs.db.eg to obtain the GO annotation and using hypergeometric test
#' (phyper) and multtest package for adjusting p-values
#' 
#' 
#' @param annotatedPeak A GRanges object or a vector of feature IDs
#' @param orgAnn Organism annotation package such as org.Hs.eg.db for human and
#' org.Mm.eg.db for mouse, org.Dm.eg.db for fly, org.Rn.eg.db for rat,
#' org.Sc.eg.db for yeast and org.Dr.eg.db for zebrafish
#' @param feature_id_type The feature type in annotatedPeak such as
#' ensembl_gene_id, refseq_id, gene_symbol or entrez_id
#' @param maxP The maximum p-value to be considered to be significant
#' @param minGOterm The minimum count in a genome for a GO term to be included
#' @param multiAdjMethod The multiple testing procedures, for details, see
#' mt.rawp2adjp in multtest package
#' @param condense Condense the results or not.
#' @param removeAncestorByPval Remove ancestor by p-value. P-value is
#' calculated by fisher exact test. If gene number in all of the children is
#' significant greater than it in parent term, the parent term will be removed
#' from the list.
#' @param keepByLevel If the shortest path from the go term to 'all' is greater
#' than the given level, the term will be removed.
#' @param subGroupComparison A logical vector to split the peaks into two 
#' groups. The enrichment analysis will compare the over-present GO terms
#' in TRUE group and FALSE group separately. The analysis will split into 
#' two steps: 1. enrichment analysis for TRUE group by hypergeometric
#' test; 2. enrichment analysis for TRUE over FALSE group
#' by Fisher's Exact test for the enriched GO terms.
#' To keep the output same format, if you want to compare FALSE vs TRUE,
#' please repeat the analysis by inverting the parameter.
#' Default is NULL.
#' @return A list with 3 elements \item{list("bp")}{ enriched biological
#' process with the following 9 variables
#' 
#' go.id:GO biological process id
#' 
#' go.term:GO biological process term
#' 
#' go.Definition:GO biological process description
#' 
#' Ontology: Ontology branch, i.e. BP for biological process
#' 
#' count.InDataset: count of this GO term in this dataset
#' 
#' count.InGenome: count of this GO term in the genome
#' 
#' pvalue: pvalue from the hypergeometric test
#' 
#' totaltermInDataset: count of all GO terms in this dataset
#' 
#' totaltermInGenome: count of all GO terms in the genome
#' 
#' } \item{list("mf")}{enriched molecular function with the following 9
#' variables
#' 
#' go.id:GO molecular function id
#' 
#' go.term:GO molecular function term
#' 
#' go.Definition:GO molecular function description
#' 
#' Ontology: Ontology branch, i.e. MF for molecular function
#' 
#' count.InDataset: count of this GO term in this dataset
#' 
#' count.InGenome: count of this GO term in the genome
#' 
#' pvalue: pvalue from the hypergeometric test
#' 
#' totaltermInDataset: count of all GO terms in this dataset
#' 
#' totaltermInGenome: count of all GO terms in the genome
#' 
#' } \item{list("cc")}{enriched cellular component the following 9 variables
#' 
#' go.id:GO cellular component id
#' 
#' go.term:GO cellular component term
#' 
#' go.Definition:GO cellular component description
#' 
#' Ontology: Ontology type, i.e. CC for cellular component
#' 
#' count.InDataset: count of this GO term in this dataset
#' 
#' count.InGenome: count of this GO term in the genome
#' 
#' pvalue: pvalue from the hypergeometric test
#' 
#' totaltermInDataset: count of all GO terms in this dataset
#' 
#' totaltermInGenome: count of all GO terms in the genome
#' 
#' }
#' @author Lihua Julie Zhu. Jianhong Ou for subGroupComparison
#' @seealso phyper, hyperGtest
#' @references Johnson, N. L., Kotz, S., and Kemp, A. W. (1992) Univariate
#' Discrete Distributions, Second Edition. New York: Wiley
#' @keywords misc
#' @export
#' @importFrom AnnotationDbi mappedkeys
#' @importFrom multtest mt.rawp2adjp
#' @importFrom AnnotationDbi Definition Ontology Term
#' @importFrom stats phyper fisher.test
#' @examples
#' 
#'   data(enrichedGO)
#'   enrichedGO$mf[1:10,]
#'   enrichedGO$bp[1:10,]
#'   enrichedGO$cc
#'   if (interactive()) {
#'      data(annotatedPeak)
#'      library(org.Hs.eg.db)
#'      library(GO.db)
#'      enriched.GO = getEnrichedGO(annotatedPeak[1:6,], 
#'                                  orgAnn="org.Hs.eg.db", 
#'                                  maxP=0.01,
#'                                  minGOterm=10,
#'                                  multiAdjMethod= NULL)
#'      dim(enriched.GO$mf)
#'      colnames(enriched.GO$mf)
#'      dim(enriched.GO$bp)
#'      enriched.GO$cc
#' }
#' 
getEnrichedGO <- function(annotatedPeak, orgAnn, 
                          feature_id_type="ensembl_gene_id", 
                          maxP=0.01,
                          minGOterm=10, multiAdjMethod=NULL,
                          condense=FALSE,
                          removeAncestorByPval=NULL,
                          keepByLevel=NULL,
                          subGroupComparison=NULL){
  if(is(annotatedPeak, "GRangesList")){
    if(length(subGroupComparison)){
      stop("subGroupComparison parameter is not supported for list of peaks.")
    }
    args <- as.list(match.call())
    res <- lapply(annotatedPeak, function(.ele){
      args$annotatedPeak <- .ele
      do.call(getEnrichedGO, args = args)
    })
  }
  stopifnot("The 'GO.db' package is required"=
              requireNamespace("GO.db", quietly = TRUE)) 
    if (missing(annotatedPeak))
    {
        stop("Missing required argument annotatedPeak!")    
    }
    if(length(multiAdjMethod)>0){
        multiAdjMethod <- match.arg(multiAdjMethod, 
                                    c("Bonferroni", "Holm", "Hochberg", 
                                      "SidakSS", "SidakSD", "BH", "BY",
                                      "ABH","TSBH"))
    }
    if (missing(orgAnn))
    {
        message("No valid organism specific GO gene mapping package 
                as argument orgAnn is passed in!")
        stop("Please refer 
             http://www.bioconductor.org/packages/release/data/annotation/ 
             for available org.xx.eg.db packages")
    }
    GOgenome = sub(".db","",orgAnn)
    if (nchar(GOgenome) <1)
    {
        message("No valid organism specific GO gene mapping package as 
                parameter orgAnn is passed in!")
        stop("Please refer 
             http://www.bioconductor.org/packages/release/data/annotation/ 
             for available org.xx.eg.db packages")
    }
    groupFALSE <- NULL
    feature_ids_FALSE <- NULL
    entrezIDs_FALSE <- NULL
    if(length(subGroupComparison)){
      if(length(subGroupComparison)!=length(annotatedPeak)){
        stop("Length of subGroupComparison should keep same as ",
             "length of annotatedPeak or check your input annotatedPeak type.")
      }
      stopifnot(is(subGroupComparison, "logical"))
      groupFALSE <- annotatedPeak[!subGroupComparison]
      annotatedPeak <- annotatedPeak[subGroupComparison]
    }
    if (inherits(annotatedPeak, what=c("GRanges")))
    {
        feature_ids = unique(as.character(annotatedPeak$feature))
        if(length(groupFALSE)){
          feature_ids_FALSE = unique(as.character(groupFALSE$feature))
        }
    }else if(is.character(annotatedPeak))
    {
        feature_ids = unique(annotatedPeak)
        if(length(groupFALSE)){
          feature_ids_FALSE = unique(groupFALSE)
        }
    }else{
        stop("annotatedPeak needs to be GRanges type with feature 
             variable holding the feature id or a character vector 
             holding the IDs of the features used to annotate the peaks!")
    }
    if (feature_id_type == "entrez_id"){
        entrezIDs <- feature_ids
        if(length(groupFALSE)){
          entrezIDs_FALSE <- feature_ids_FALSE
        }
    }else{
        cov2EntrezID <- function(IDs, orgAnn, ID_type="ensembl_gene_id"){
            GOgenome = sub(".db","",orgAnn)
            orgAnn <- switch(ID_type,
                             ensembl_gene_id="ENSEMBL2EG",
                             gene_symbol="SYMBOL2EG",
                             refseq_id="REFSEQ2EG",
                             "BAD_NAME")
            if(orgAnn=="BAD_NAME")
                stop("Currently only the following type of IDs are supported:
                     ensembl_gene_id, refseq_id and gene_symbol!")
            orgAnn <- get(paste(GOgenome, orgAnn, sep=""))
            if(!is(orgAnn, "AnnDbBimap")){
                stop("orgAnn is not a valid annotation dataset! 
                     For example, orgs.Hs.eg.db package for human and 
                     the org.Mm.eg.db package for mouse.")
            }
            IDs <- unique(IDs[!is.na(IDs)])
            ids <- mget(IDs, orgAnn, ifnotfound=NA)
            ids <- lapply(ids, `[`, 1)
            ids <- unlist(ids)
            ids <- unique(ids[!is.na(ids)])
            ids
        }
        entrezIDs <- cov2EntrezID(feature_ids, orgAnn, feature_id_type)
        if(length(groupFALSE)){
          entrezIDs_FALSE <- cov2EntrezID(feature_ids_FALSE,
                                          orgAnn, feature_id_type)
        }
    }
    if(length(entrezIDs)<2){
        stop("The number of gene is less than 2. 
             Please double check your feature_id_type.")
    }
    goAnn <- get(paste(GOgenome,"GO", sep=""))
    mapped_genes <- mappedkeys(goAnn)
    totalN.genes=length(unique(mapped_genes))
    thisN.genes = length(unique(entrezIDs))
    #xx <- as.list(goAnn[mapped_genes])
    
    xx <- mget(mapped_genes, goAnn, ifnotfound=NA)
    all.GO <- cbind(matrix(unlist(unlist(xx)),ncol=3,byrow=TRUE),
                    rep(names(xx), elementNROWS(xx)))
    all.GO <- unique(all.GO)## incalse the database is not unique
    this.GO <- all.GO[all.GO[, 4] %in% entrezIDs, , drop=FALSE]
    FALSE.GO <- all.GO[all.GO[, 4] %in% entrezIDs_FALSE, , drop=FALSE]
    
    addAnc <- function(go.ids,
                       onto=c("bp", "cc", "mf")){
        ##replace the function addAncestors
        GOIDs <- unique(as.character(go.ids[, 1]))
        empty <- matrix(nrow=0, ncol=2)
        colnames(empty) <- c("go.id", "EntrezID")
        if(length(GOIDs)<1){
            return(empty)
        }
        onto <- match.arg(onto)
        Ancestors <- switch(onto, 
                            bp=mget(GOIDs, GO.db::GOBPANCESTOR, ifnotfound = NA),
                            cc=mget(GOIDs, GO.db::GOCCANCESTOR, ifnotfound = NA),
                            mf=mget(GOIDs, GO.db::GOMFANCESTOR, ifnotfound = NA))
        Ancestors <- unique(unlist(Ancestors))
        Ancestors <- Ancestors[Ancestors!="all" & !is.na(Ancestors)]
        if(length(Ancestors)>0){
            children <- go.ids[,c(1, 4), drop=FALSE]
            children <- unique(children)
            children.s <- split(children[, 2], children[, 1])
            Ancestors <- switch(onto,
                               bp=mget(Ancestors, 
                                       GO.db::GOBPOFFSPRING, 
                                       ifnotfound = NA),
                               cc=mget(Ancestors, 
                                       GO.db::GOCCOFFSPRING, 
                                       ifnotfound = NA),
                               mf=mget(Ancestors, 
                                       GO.db::GOMFOFFSPRING, 
                                       ifnotfound = NA))
            Ancestors <- cbind(Ancestor=rep(names(Ancestors),
                                           elementNROWS(Ancestors)),
                              child=unlist(Ancestors))
            Ancestors <- 
                Ancestors[Ancestors[, "child"] %in% names(children.s), 
                          , drop=FALSE]
            temp <- cbind(Ancestor=children[, 1], child=children[, 1])
            Ancestors <- rbind(Ancestors, temp)
            Ancestors <- unique(Ancestors)
            Ancestors.ezid <- children.s[Ancestors[, "child"]]
            temp <- cbind(rep(Ancestors[, "Ancestor"], 
                                    elementNROWS(Ancestors.ezid)),
                          unlist(Ancestors.ezid))
            temp <- unique(temp) #time....
            if (length(temp) <3){
                re <- unique(children)
            }else{
                re <- temp[!is.na(temp[,1]) & temp[,1] != "",]
            }
        }else{
            re <- unique(cbind(as.character(go.ids[,1]), go.ids[,4]))
        }
        colnames(re) <- c("go.id", "EntrezID")
        re
    }
    
    bp.go.this.withEntrez = addAnc(this.GO[this.GO[,3]=="BP",],"bp")
    cc.go.this.withEntrez = addAnc(this.GO[this.GO[,3]=="CC",], "cc")
    mf.go.this.withEntrez = addAnc(this.GO[this.GO[,3]=="MF",],"mf")
    
    bp.go.all.withEntrez  = addAnc(all.GO[all.GO[,3]=="BP",], "bp")
    cc.go.all.withEntrez = addAnc(all.GO[all.GO[,3]=="CC",], "cc")
    mf.go.all.withEntrez = addAnc(all.GO[all.GO[,3]=="MF",], "mf")
    
    bp.go.this  = bp.go.this.withEntrez[,1]
    cc.go.this = cc.go.this.withEntrez[,1]
    mf.go.this = mf.go.this.withEntrez[,1]
    
    bp.go.all  = bp.go.all.withEntrez[,1]
    cc.go.all = cc.go.all.withEntrez[,1]
    mf.go.all = mf.go.all.withEntrez[,1]
    
    total.mf = length(mf.go.all)
    total.cc = length(cc.go.all)
    total.bp = length(bp.go.all)
    this.mf = length(mf.go.this)
    this.cc = length(cc.go.this)
    this.bp = length(bp.go.this)
    
    this.bp.count <- table(as.character(bp.go.this[bp.go.this!=""]))
    this.mf.count <- table(as.character(mf.go.this[mf.go.this!=""]))
    this.cc.count <- table(as.character(cc.go.this[cc.go.this!=""]))
    
    all.bp.count <- table(as.character(bp.go.all[bp.go.all!=""]))
    all.mf.count <- table(as.character(mf.go.all[mf.go.all!=""]))
    all.cc.count <- table(as.character(cc.go.all[cc.go.all!=""]))
    
    if(length(groupFALSE)){
      bp.go.FALSE.withEntrez = addAnc(FALSE.GO[FALSE.GO[,3]=="BP",],"bp")
      cc.go.FALSE.withEntrez = addAnc(FALSE.GO[FALSE.GO[,3]=="CC",], "cc")
      mf.go.FALSE.withEntrez = addAnc(FALSE.GO[FALSE.GO[,3]=="MF",],"mf")
      
      bp.go.FALSE  = bp.go.FALSE.withEntrez[,1]
      cc.go.FALSE = cc.go.FALSE.withEntrez[,1]
      mf.go.FALSE = mf.go.FALSE.withEntrez[,1]
      
      FALSE.mf = length(mf.go.FALSE)
      FALSE.cc = length(cc.go.FALSE)
      FALSE.bp = length(bp.go.FALSE)
      
      FALSE.bp.count <- table(as.character(bp.go.FALSE[bp.go.FALSE!=""]))
      FALSE.mf.count <- table(as.character(mf.go.FALSE[mf.go.FALSE!=""]))
      FALSE.cc.count <- table(as.character(cc.go.FALSE[cc.go.FALSE!=""]))
    }
    
    hyperGT <- function(alltermcount, thistermcount, 
                        totaltermInGenome, totaltermInPeakList){
        m <- as.numeric(alltermcount[names(thistermcount)])
        q <- as.numeric(thistermcount)
        n <- as.numeric(totaltermInGenome)
        k <- as.numeric(totaltermInPeakList)
        pvalue <- phyper(q-1, m, n-m, k, lower.tail = FALSE, log.p = FALSE)
        data.frame(go.id=names(thistermcount), 
                   count.InDataset=q,
                   count.InGenome=m, 
                   pvalue=pvalue,
                   totaltermInDataset=k, 
                   totaltermInGenome=n)
    }
    
    bp.selected = hyperGT(all.bp.count,
                             this.bp.count, 
                             total.bp, 
                             this.bp)
    mf.selected = hyperGT(all.mf.count,
                             this.mf.count, 
                             total.mf, 
                             this.mf)
    cc.selected = hyperGT(all.cc.count,
                             this.cc.count, 
                             total.cc, 
                             this.cc)
    
    if (length(multiAdjMethod)<1){
        bp.s = 
            bp.selected[
                as.numeric(as.character(bp.selected[,4]))<maxP & 
                    as.numeric(as.character(bp.selected[,3]))>=minGOterm,]
        mf.s = 
            mf.selected[
                as.numeric(as.character(mf.selected[,4]))<maxP & 
                    as.numeric(as.character(mf.selected[,3]))>=minGOterm,]
        cc.s = 
            cc.selected[
                as.numeric(as.character(cc.selected[,4]))<maxP & 
                    as.numeric(as.character(cc.selected[,3]))>=minGOterm,]
    }else{
        procs = c(multiAdjMethod)
        res <- mt.rawp2adjp(as.numeric(as.character(bp.selected[,4])), procs)
        adjp = unique(res$adjp)
        colnames(adjp)[1] = colnames(bp.selected)[4]
        colnames(adjp)[2] = paste(multiAdjMethod, "adjusted.p.value", sep=".")
        bp.selected[,4] = as.numeric(as.character(bp.selected[,4]))
        bp1 = merge(bp.selected, adjp, all.x=TRUE)
        
        res <- mt.rawp2adjp(as.numeric(as.character(mf.selected[,4])), procs)
        adjp = unique(res$adjp)
        colnames(adjp)[1] = colnames(mf.selected)[4]
        colnames(adjp)[2] = paste(multiAdjMethod, "adjusted.p.value", sep=".")
        mf.selected[,4] = as.numeric(as.character(mf.selected[,4]))
        mf1 = merge(mf.selected, adjp, all.x=TRUE)
        
        res <- mt.rawp2adjp(as.numeric(as.character(cc.selected[,4])), procs)
        adjp = unique(res$adjp)
        colnames(adjp)[1] = colnames(cc.selected)[4]
        colnames(adjp)[2] = paste(multiAdjMethod, "adjusted.p.value", sep=".")
        cc.selected[,4] = as.numeric(as.character(cc.selected[,4]))
        cc1 = merge(cc.selected, adjp, all.x=TRUE)
        
        bp.s = bp1[as.numeric(as.character(bp1[,dim(bp1)[2]]))<maxP &  
                       !is.na(bp1[,dim(bp1)[2]]) & 
                       as.numeric(as.character(bp1[,4]))>=minGOterm,]
        mf.s = mf1[as.numeric(as.character(mf1[,dim(mf1)[2]]))<maxP & 
                       !is.na(mf1[,dim(mf1)[2]]) & 
                       as.numeric(as.character(mf1[,4]))>=minGOterm,]
        cc.s = cc1[as.numeric(as.character(cc1[,dim(cc1)[2]]))<maxP & 
                       !is.na(cc1[,dim(cc1)[2]]) & 
                       as.numeric(as.character(cc1[,4]))>=minGOterm,]
    }
    
    annoTerms <- function(goids){
        if(length(goids)<1){
            goterm <- matrix(ncol=4)
        }else{
            goids <- as.character(goids)
            terms <- Term(goids)
            definition <- Definition(goids)
            ontology <- Ontology(goids)
            goterm <- cbind(goids, 
                            terms[match(goids, names(terms))],
                            definition[match(goids, names(definition))],
                            ontology[match(goids, names(ontology))])
        }
        colnames(goterm) <- c("go.id", "go.term", "Definition", "Ontology")
        rownames(goterm) <- NULL
        goterm
    }
    goterm.bp <- annoTerms(bp.s$go.id)
    goterm.mf <- annoTerms(mf.s$go.id)
    goterm.cc <- annoTerms(cc.s$go.id)
    
    bp.selected1 = merge(goterm.bp, bp.s,by="go.id")
    mf.selected1 = merge(goterm.mf, mf.s,by="go.id")
    cc.selected1 = merge(goterm.cc, cc.s,by="go.id")
    
    if(length(groupFALSE)){ ## add counts for FALSE group and do Fisher's exact test
      bp.selected1$count.InBackgroundDataset = FALSE.bp.count[bp.selected1$go.id]
      mf.selected1$count.InBackgroundDataset = FALSE.mf.count[mf.selected1$go.id]
      cc.selected1$count.InBackgroundDataset = FALSE.cc.count[cc.selected1$go.id]
      
      multiFisher <- function(.data){
        rowFisher <- function(x, ...){
          return(fisher.test(matrix(x, nrow = 2, byrow = TRUE), ...)$p.value)
        }
        .data0 <- .data[, c("count.InDataset", "count.InBackgroundDataset")]
        .data0 <- as.matrix(.data0)
        .data0[is.na(.data0)] <- 0
        .data0 <- cbind(.data0, .data[, "count.InGenome"] - .data0)
        apply(.data0, 1, rowFisher, alternative = "greater")
      }
      bp.selected1$dataset.vs.background.pval <- multiFisher(bp.selected1)
      mf.selected1$dataset.vs.background.pval <- multiFisher(mf.selected1)
      cc.selected1$dataset.vs.background.pval <- multiFisher(cc.selected1)
      
      bp.selected1 <-
        bp.selected1[bp.selected1$dataset.vs.background.pval<maxP, , drop=FALSE]
      mf.selected1 <-
        mf.selected1[mf.selected1$dataset.vs.background.pval<maxP, , drop=FALSE]
      cc.selected1 <-
        cc.selected1[cc.selected1$dataset.vs.background.pval<maxP, , drop=FALSE]
    }
    
    if(condense){
        bp.go.this.withEntrez <- 
            condenseMatrixByColnames(bp.go.this.withEntrez, iname="go.id")
        mf.go.this.withEntrez <- 
            condenseMatrixByColnames(mf.go.this.withEntrez, iname="go.id")
        cc.go.this.withEntrez <- 
            condenseMatrixByColnames(cc.go.this.withEntrez, iname="go.id")
    }
    
    bp.selected = merge(bp.selected1, bp.go.this.withEntrez)
    mf.selected = merge(mf.selected1, mf.go.this.withEntrez)
    cc.selected = merge(cc.selected1, cc.go.this.withEntrez)
    
    if(length(removeAncestorByPval)>0){
      if(is.numeric(removeAncestorByPval[1]) || 
         is.integer(removeAncestorByPval[1])){
        bp.selected <- removeAncestor(bp.selected, onto="BP", 
                                      cutoffPvalue=removeAncestorByPval)
        mf.selected <- removeAncestor(mf.selected, onto="MF", 
                                      cutoffPvalue=removeAncestorByPval)
        cc.selected <- removeAncestor(cc.selected, onto="CC", 
                                      cutoffPvalue=removeAncestorByPval)
      }
    }
    
    if(length(keepByLevel)>0){
      if(is.numeric(keepByLevel[1]) || is.integer(keepByLevel[1])){
        bp.selected <- filterByLevel(bp.selected, onto="BP", level=keepByLevel)
        mf.selected <- filterByLevel(mf.selected, onto="MF", level=keepByLevel)
        cc.selected <- filterByLevel(cc.selected, onto="CC", level=keepByLevel)
      }
    }
    list(bp=bp.selected, mf=mf.selected, cc=cc.selected)
}
