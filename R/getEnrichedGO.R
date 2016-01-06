getEnrichedGO <- function(annotatedPeak, orgAnn, 
                          feature_id_type="ensembl_gene_id", 
                          maxP=0.01, multiAdj=FALSE, 
                          minGOterm=10, multiAdjMethod="",
                          condense=FALSE){    
    if (missing(annotatedPeak))
    {
        stop("Missing required argument annotatedPeak!")    
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
    if (inherits(annotatedPeak, what=c("RangedData", "GRanges")))
    {
        feature_ids = unique(as.character(annotatedPeak$feature))
    }
    else if (class(annotatedPeak)  ==  "character")
    {
        feature_ids = unique(annotatedPeak)
    }
    else
    {
        stop("annotatedPeak needs to be GRanges type with feature 
             variable holding the feature id or a character vector 
             holding the IDs of the features used to annotate the peaks!")
    }
    if (feature_id_type == "entrez_id")
    {
        entrezIDs <- feature_ids
    }
    else
    {
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
            if(class(orgAnn) != "AnnDbBimap"){
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
    }
    if(length(entrezIDs)<2){
        stop("The number of gene is less than 2. Please double check your feature_id_type.")
    }
    
    goAnn <- get(paste(GOgenome,"GO", sep=""))
    mapped_genes <- mappedkeys(goAnn)
    totalN.genes=length(unique(mapped_genes))
    thisN.genes = length(unique(entrezIDs))
    #xx <- as.list(goAnn[mapped_genes])
    
    xx <- mget(mapped_genes, goAnn, ifnotfound=NA)
    all.GO <- cbind(matrix(unlist(unlist(xx)),ncol=3,byrow=TRUE),
                    rep(names(xx), elementLengths(xx)))
    
    this.GO <- all.GO[all.GO[, 4] %in% entrezIDs, , drop=FALSE]
    
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
                            bp=mget(GOIDs, GOBPANCESTOR, ifnotfound = NA),
                            cc=mget(GOIDs, GOCCANCESTOR, ifnotfound = NA),
                            mf=mget(GOIDs, GOMFANCESTOR, ifnotfound = NA))
        Ancestors <- unique(unlist(Ancestors))
        Ancestors <- Ancestors[Ancestors!="all" & !is.na(Ancestors)]
        if(length(Ancestors)>0){
            children <- go.ids[,c(1, 4), drop=FALSE]
            children <- unique(children)
            children.s <- split(children[, 2], children[, 1])
            Ancestors <- switch(onto,
                               bp=mget(Ancestors, 
                                       GOBPOFFSPRING, 
                                       ifnotfound = NA),
                               cc=mget(Ancestors, 
                                       GOCCOFFSPRING, 
                                       ifnotfound = NA),
                               mf=mget(Ancestors, 
                                       GOMFOFFSPRING, 
                                       ifnotfound = NA))
            Ancestors <- cbind(Ancestor=rep(names(Ancestors),
                                           elementLengths(Ancestors)),
                              child=unlist(Ancestors))
            Ancestors <- 
                Ancestors[Ancestors[, "child"] %in% names(children.s), 
                          , drop=FALSE]
            temp <- cbind(Ancestor=children[, 1], child=children[, 1])
            Ancestors <- rbind(Ancestors, temp)
            Ancestors <- unique(Ancestors)
            Ancestors.ezid <- children.s[Ancestors[, "child"]]
            temp <- cbind(rep(Ancestors[, "Ancestor"], 
                                    elementLengths(Ancestors.ezid)),
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
    
    if (multiAdj == FALSE | multiAdjMethod == "")
    {
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
    }
    else
    {
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
    
    list(bp=bp.selected, mf=mf.selected, cc=cc.selected)
}