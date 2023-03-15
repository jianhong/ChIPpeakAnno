#' Obtain enriched PATH that near the peaks
#'
#' Obtain enriched PATH that are near the peaks using path package such as
#' reactome.db and path mapping package such as org.Hs.db.eg to obtain the path
#' annotation and using hypergeometric test (phyper) and multtest package for
#' adjusting p-values
#'
#'
#' @param annotatedPeak GRanges such as data(annotatedPeak) or a vector of
#' feature IDs
#' @param orgAnn organism annotation package such as org.Hs.eg.db for human and
#' org.Mm.eg.db for mouse, org.Dm.eg.db for fly, org.Rn.eg.db for rat,
#' org.Sc.eg.db for yeast and org.Dr.eg.db for zebrafish
#' @param pathAnn pathway annotation package such as KEGG.db (deprecated), reactome.db, KEGGREST
#' @param feature_id_type the feature type in annotatedPeakRanges such as
#' ensembl_gene_id, refseq_id, gene_symbol or entrez_id
#' @param maxP maximum p-value to be considered to be significant
#' @param minPATHterm minimum count in a genome for a path to be included
#' @param multiAdjMethod multiple testing procedures, for details, see
#' mt.rawp2adjp in multtest package
#' @param subGroupComparison A logical vector to split the peaks into two
#' groups. The enrichment analysis will compare the over-present GO terms
#' in TRUE group and FALSE group separately. The analysis will split into
#' two steps: 1. enrichment analysis for TRUE group by hypergeometric
#' test; 2. enrichment analysis for TRUE over FALSE group
#' by Fisher's Exact test for the enriched GO terms.
#' To keep the output same format, if you want to compare FALSE vs TRUE,
#' please repeat the analysis by inverting the parameter.
#' Default is NULL.
#' @return A dataframe of enriched path with the following variables.
#' \item{path.id}{KEGG PATH ID} \item{EntrezID}{Entrez ID}
#' \item{count.InDataset}{count of this PATH in this dataset}
#' \item{count.InGenome}{count of this PATH in the genome} \item{pvalue}{pvalue
#' from the hypergeometric test} \item{totaltermInDataset}{count of all PATH in
#' this dataset} \item{totaltermInGenome}{count of all PATH in the genome}
#' \item{PATH}{PATH name}
#' @author Jianhong Ou, Kai Hu
#' @seealso phyper, hyperGtest
#' @references Johnson, N. L., Kotz, S., and Kemp, A. W. (1992) Univariate
#' Discrete Distributions, Second Edition. New York: Wiley
#' @keywords misc
#' @export
#' @importFrom AnnotationDbi mappedkeys
#' @importFrom multtest mt.rawp2adjp
#' @importFrom KEGGREST keggGet keggLink
#' @examples
#'
#' if (interactive()||Sys.getenv("USER")=="jianhongou") {
#' data(annotatedPeak)
#' library(org.Hs.eg.db)
#' library(reactome.db)
#' enriched.PATH = getEnrichedPATH(annotatedPeak, orgAnn="org.Hs.eg.db",
#'                  feature_id_type="ensembl_gene_id",
#'                  pathAnn="reactome.db", maxP=0.01,
#'                  minPATHterm=10, multiAdjMethod=NULL)
#'  head(enriched.PATH)
#'  enrichedKEGG = getEnrichedPATH(annotatedPeak, orgAnn="org.Hs.eg.db",
#'                  feature_id_type="ensembl_gene_id",
#'                  pathAnn="KEGGREST", maxP=0.01,
#'                  minPATHterm=10, multiAdjMethod=NULL)
#'  enrichmentPlot(enrichedKEGG)
#' }
#'
getEnrichedPATH <- function(annotatedPeak, orgAnn, pathAnn,
                            feature_id_type="ensembl_gene_id",
                            maxP=0.01, minPATHterm=10, multiAdjMethod=NULL,
                            subGroupComparison=NULL)
{
    if (missing(annotatedPeak)){
        stop("Missing required argument annotatedPeak!")
    }
    # if (missing(feature_id_type)) {
    #     stop("Missing required argument feature_id_type!")
    # }
    if(length(multiAdjMethod)>0){
        multiAdjMethod <- match.arg(multiAdjMethod,
                                    c("Bonferroni", "Holm", "Hochberg",
                                      "SidakSS", "SidakSD", "BH", "BY",
                                      "ABH","TSBH"))
    }
    if (!grepl("^org\\...\\.eg\\.db",orgAnn)){
        message("No valid organism specific PATH gene mapping package as
                parameter orgAnn is passed in!")
        stop("Please refer
             http://www.bioconductor.org/packages/release/data/annotation/
             for available org.xx.eg.db packages")
    }
    if(!isNamespaceLoaded(orgAnn)){
        stop(paste("Need to load", orgAnn,
                   "before using getEnrichedPATH. Try \n\"library(",
                   orgAnn,")\""))
    }
    if (missing(pathAnn)){
        stop("Missing required argument pathAnn. \n
             pathAnn is the database with annotation object that maps
Entrez Gene to pathway identifies named as xxxxxEXTID2PATHID
             and pathway identifies to pathway names named as
             xxxxxPATHID2NAME.")
    }
    if(!isNamespaceLoaded(pathAnn)){
        stop(paste("Need to load",pathAnn,
                   "before using getEnrichedPATH. Try \n\"library(",
                   pathAnn,")\""))
    }

    groupFALSE <- NULL
    feature_ids_FALSE <- NULL
    entrezIDs_FALSE <- NULL
    if(length(subGroupComparison)){
        if(length(subGroupComparison)!=length(annotatedPeak)){
            stop("Length of subGroupComparison should keep same as",
                 "length of annotatedPeak or check your input annotatedPeak type.")
        }
        stopifnot(is(subGroupComparison, "logical"))
        groupFALSE <- annotatedPeak[!subGroupComparison]
        annotatedPeak <- annotatedPeak[subGroupComparison]
    }
    if (inherits(annotatedPeak, what=c("GRanges"))){
        feature_ids = unique(as.character(annotatedPeak$feature))
        if(length(groupFALSE)){
            feature_ids_FALSE = unique(as.character(groupFALSE$feature))
        }
    }else if (is.character(annotatedPeak)){
        feature_ids = unique(annotatedPeak)
        if(length(groupFALSE)){
            feature_ids_FALSE = unique(groupFALSE)
        }
    }else{
        stop("annotatedPeak needs to be GRanges with feature variable
             holding the feature id or a character vector holding the IDs of
             the features used to annotate the peaks!")
    }
    if (feature_id_type == "entrez_id"){
        entrezIDs <- feature_ids
        if(length(groupFALSE)){
            entrezIDs_FALSE <- feature_ids_FALSE
        }
    }else{
        entrezIDs <- convert2EntrezID(feature_ids, orgAnn, feature_id_type)
        if(length(groupFALSE)){
            entrezIDs_FALSE <- convert2EntrezID(feature_ids_FALSE,
                                            orgAnn, feature_id_type)
        }
    }
    if(length(entrezIDs)<2){
        stop('The number of gene is less than 2. Please double check your feature_id_type. If using "ensembl_gene_id", do not include the version number. E.g. "ENSG00000139618" is okay while "ENSG00000139618.2" is not compatible. Use `ensembl_gene_id <- sub("\\.[^.]+$", "", ensembl_gene_id)` to batch remove the version number if needed.')
    }
    if (pathAnn %in% c("KEGG.db")) {
      stop('`pathAnn = "KEGG.db"` has been deprecated, use `pathAnn = "KEGGREST"` instead.')
    }
    if (pathAnn %in% c("reactome.db", "KEGG.db")) {
        extid2path<- paste(gsub(".db$","",pathAnn),"EXTID2PATHID", sep="")
        path2name<- paste(gsub(".db$","",pathAnn),"PATHID2NAME", sep="")
        if(length(objects(paste("package",pathAnn,sep=":"),
                          pattern=extid2path))!=1 &
           length(objects(paste("package",pathAnn,sep=":"),
                          pattern=path2name))!=1 ){
            stop("argument pathAnn is not the annotation data with objects named
             as xxxxxEXTID2PATHID and/or xxxxxPATHID2NAME")
        }
        extid2path <- get(extid2path)
        mapped_genes <- mappedkeys(extid2path)
        #get all the entrez_ids in the species
        mapped_genes <-
            mapped_genes[mapped_genes %in%
                             mappedkeys(get(paste(gsub(".db","",orgAnn),
                                                  "SYMBOL",sep="")))]
        totalN.genes=length(unique(mapped_genes))
        thisN.genes = length(unique(entrezIDs))
        xx <- as.list(extid2path[mapped_genes])
    } else if (pathAnn == "KEGGREST") {
    # for KEGGREST db
        organismKEGGREST <- .findKEGGRESTOrganismName(orgAnn)
        EGID2PATHID <- keggLink("pathway", organismKEGGREST)
        # get rid of the leading organismKEGGREST in front of the EntrezID and the "path" in front of the PATHID
        names(EGID2PATHID) <- unlist(lapply(strsplit(names(EGID2PATHID), ":"), "[", 2))
        EGID2PATHID <- unlist(lapply(strsplit(EGID2PATHID, ":"), "[", 2))

        mapped_genes2 <- names(keggLink("pathway", organismKEGGREST))
        mapped_genes2 <- unique(unlist(lapply(strsplit(mapped_genes2, ":"), "[", 2)))

        mapped_genes <- mapped_genes2[mapped_genes2 %in%
                            mappedkeys(get(paste(gsub(".db","",orgAnn), "SYMBOL", sep="")))]

        totalN.genes=length(unique(mapped_genes))
        thisN.genes = length(unique(entrezIDs))

        xx <- with(stack(as.list(EGID2PATHID)), split(values, ind))
    }

    all.PATH <- do.call(rbind, lapply(mapped_genes,function(x1)
    {
        temp = unlist(xx[names(xx) ==x1])
        if (length(temp) >0)
        {
            temp1 =matrix(temp,ncol=1,byrow=TRUE)
            cbind(temp1,rep(x1,dim(temp1)[1]))
        }
    }))
    this.PATH <- do.call(rbind, lapply(entrezIDs,function(x1)
    {
        temp = unlist(xx[names(xx) ==x1])
        if (length(temp) >0)
        {
            temp1 =matrix(temp,ncol=1,byrow=TRUE)
            cbind(temp1,rep(x1,dim(temp1)[1]))
        }
    }))


    all.PATH <- unique(all.PATH)## in case the database is not unique
    this.PATH <- unique(this.PATH)
    if (is.null(all.PATH) | is.null(this.PATH)) {
    	return("No enriched pathway found!")
    }
    colnames(all.PATH)<-c("path.id","entrez_id")
    colnames(this.PATH)<-c("path.id","entrez_id")
    path.all<-as.character(all.PATH[,"path.id"])
    path.this<-as.character(this.PATH[,"path.id"])

    total = length(path.all)
    this = length(path.this)

    all.count = getUniqueGOidCount(as.character(path.all[path.all!=""]))
    this.count = getUniqueGOidCount(as.character(path.this[path.this!=""]))

    if(length(groupFALSE)){
        FALSE.PATH <- do.call(rbind, lapply(entrezIDs_FALSE,function(x1)
        {
            temp = unlist(xx[names(xx) ==x1])
            if (length(temp) >0)
            {
                temp1 =matrix(temp,ncol=1,byrow=TRUE)
                cbind(temp1,rep(x1,dim(temp1)[1]))
            }
        }))
        FALSE.PATH <- unique(FALSE.PATH)
        colnames(FALSE.PATH)<-c("path.id","entrez_id")
        path.FALSE<-as.character(FALSE.PATH[,"path.id"])
        FALSE.count = getUniqueGOidCount(as.character(path.FALSE[path.FALSE!=""]))
        names(FALSE.count[[2]]) <- FALSE.count[[1]]
        FALSE.count <- FALSE.count[[2]]
    }

    selected = hyperGtest(all.count,this.count, total, this)

    selected = data.frame(selected)

    colnames(selected) = c("path.id", "count.InDataset", "count.InGenome",
                           "pvalue", "totaltermInDataset", "totaltermInGenome")
		
    annoTerms <- function(termids){
        if(length(termids)<1){
            goterm <- matrix(ncol=2)
        }else{
            if (pathAnn %in% c("reactome.db")) {
                termids <- as.character(termids)
                terms <- xget(termids, get(sub(".db", "PATHID2NAME", pathAnn)))
                goterm <- cbind(termids,
                                terms[match(termids, names(terms))])
            } else if (pathAnn %in% c("KEGG.db", "KEGGREST")) {
                # Must get rid of the leading organismKEGGREST of path.id
                # if using KEGG.db or KEGGREST
                organismKEGGREST <- .findKEGGRESTOrganismName(orgAnn)
                termidsPreRemoved <- sub(organismKEGGREST, "", termids)

                if (pathAnn == "KEGG.db") {
                    terms <- xget(termidsPreRemoved, get(sub(".db", "PATHID2NAME", pathAnn)))
                    goterm <- cbind(termids,
                                    terms[match(termidsPreRemoved, names(terms))])
                } else if (pathAnn  == "KEGGREST") {
                    getPathName <- function(pathid) {
                        namePath <- tryCatch(
                            res <- keggGet(pathid)[[1]]$NAME,
                            error = function(c) res <- "NA"
                        )
                        names(namePath) <- pathid
                        namePath
                    }
                    terms <- unlist(lapply(termids, getPathName))
                    goterm <- cbind(termids,
                                    terms[match(termids, names(terms))])
                }
            }
        }
        colnames(goterm) <- c("path.id", "path.term")
        rownames(goterm) <- NULL
        goterm
    }
    annoterm <- annoTerms(selected$path.id)
    selected <- merge(annoterm, selected, by="path.id")

    if (is.null(multiAdjMethod)){
        s = selected[as.numeric(as.character(selected[,"pvalue"]))<maxP &
                         as.numeric(as.character(selected[, "count.InGenome"]))>=minPATHterm,]
    }else{
        procs = c(multiAdjMethod)
        res <- mt.rawp2adjp(as.numeric(as.character(selected[,"pvalue"])), procs)
        adjp = unique(res$adjp)
        colnames(adjp)[1] = "pvalue"
        colnames(adjp)[2] = paste(multiAdjMethod, "adjusted.p.value", sep=".")
        selected[,"pvalue"] = as.numeric(as.character(selected[,"pvalue"]))
        bp1 = merge(selected, adjp, all.x=TRUE)

        s = bp1[as.numeric(as.character(bp1[,dim(bp1)[2]]))<maxP &
                    !is.na(bp1[,dim(bp1)[2]]) &
                    as.numeric(as.character(bp1[,"count.InGenome"]))>=minPATHterm,]
    }

    selected = merge(this.PATH, s)

    if(length(groupFALSE)){ ## add counts for FALSE group and do Fisher's exact test
        selected$count.InBackgroundDataset = FALSE.count[selected$path.id]

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
        selected$dataset.vs.background.pval <- multiFisher(selected)

        selected <-
            selected[selected$dataset.vs.background.pval<maxP, , drop=FALSE]
    }


    selected
    }


.findKEGGRESTOrganismName <- function(name) {
    # Find the organism names used in the KEGGREST db, tested for Hs, Mm, Dm, Rn, Sc, and Dr
    #acroName <- egOrgMap(orgName)
    #KEGGOrganismList <- keggList("organism")
    #ad <- adist(acroName, KEGGOrganismList[, "species"])[1, ]
    #KEGGOrganism <- KEGGOrganismList[which.min(ad), "organism"][[1]]
    if(!is(name, "character"))
        stop("class of input organism should be character")
    organism <- c(
        "org.Ag.eg.db"="aga",
        "org.At.eg.db"="ath",
        "org.Bt.eg.db"="bta",
        "org.Ce.eg.db"="cel",
        "org.Cf.eg.db"="cfa",
        "org.Dm.eg.db"="dme",
        "org.Dr.eg.db"="dre",
        "org.EcK12.eg.db"="eco",
        "org.EcSakai.eg.db"="ecs",
        "org.Gg.eg.db"="gga",
        "org.Hs.eg.db"="hsa",
        "org.Mm.eg.db"="mmu",
        "org.Mmu.eg.db"="mcc",
        "org.Pf.plasmo.db"="pfa",
        "org.Pt.eg.db"="ptr",
        "org.Rn.eg.db"="rno",
        "org.Sc.sgd.db"="sce",
        "org.Sco.eg.db"="sco",
        "org.Ss.eg.db"="ssc",
        "org.Tgondii.eg.db"="tgo",
        "org.Xl.eg.db"="xla")
    if(name %in% names(organism)){
        return(organism[name])
    }else{
        if(name %in% organism)
            return(names(organism)[organism==name])
    }
    ##claculate the string distance, need utils package
    org <- as.character(c(names(organism), organism))
    dis <- adist(name, org)
    org <- org[dis==min(dis)]
    stop(paste("You input \"", name, "\". Do you mean \"", org, "\"?", sep=""))
}
