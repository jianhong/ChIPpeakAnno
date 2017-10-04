addGeneIDs<-function(annotatedPeak, orgAnn, IDs2Add=c("symbol"), 
                     feature_id_type="ensembl_gene_id",
                     silence=TRUE, 
                     mart){
    if (missing(annotatedPeak))
    {
        stop("Missing required argument annotatedPeak!",
             call.=FALSE)
    }
    if(missing(orgAnn) & missing(mart)){
        stop('no annotation database selected',
             call.=FALSE)
    }
    
    if(class(annotatedPeak)=="RangedData"){
        annotatedPeak <- RangedData2GRanges(annotatedPeak)
    }
    feature_ids <- switch(class(annotatedPeak),
                          GRanges=unique(annotatedPeak$feature),
                          character=unique(annotatedPeak),
                          stop("annotatedPeak needs to be GRanges type with 
                               feature variable holding the feature id or a 
                               character vector holding the IDs of the features 
                               used to annotate the peaks!",call.=FALSE))
    
    feature_ids <- feature_ids[!is.na(feature_ids)]
    feature_ids <- feature_ids[feature_ids!=""]
    if (length(feature_ids) == 0)
    {
        stop("There is no feature column in annotatedPeak or 
             annotatedPeak has size 0!",call.=FALSE)
        
    }
    if(!missing(orgAnn)){
        if(class(orgAnn)=="OrgDb"){
            orgAnn <- deparse(substitute(orgAnn))
        }
        if(class(orgAnn)!="character"){
            stop("orgAnn must be a character.")
        }
        if(!grepl(".eg.db",orgAnn,ignore.case=TRUE)){
            stop('Annotation database must be *.eg.db',call.=FALSE)
        }
        is.installed <- function(orgpkg) 
            is.element(orgpkg, installed.packages()[,1])
        if(!is.installed(orgAnn)){
            biocLite(pkgs=orgAnn, 
                     suppressUpdates=TRUE, 
                     suppressAutoUpdate=TRUE)
        }
        if(!library(orgAnn, 
                    character.only=TRUE, 
                    logical.return=TRUE)){
            if(!silence) 
                message("No valid gene mapping package as 
                        argument orgAnn is passed in!")
            stop("Please refer 
                 http://www.bioconductor.org/packages/release/data/annotation/ 
                 for available org.xx.eg.db packages")
        }
        #	require(orgAnn,character.only = TRUE)
        orgAnn<-sub("\\.db$","",orgAnn,ignore.case=TRUE)
        #get Entrez ID::entrezIDs
        if(feature_id_type=="entrez_id"){
            m_ent<-as.data.frame(feature_ids, stringsAsFactors=FALSE)
            colnames(m_ent)<-c("entrez_id")
        }else{
            prefix<-switch(feature_id_type,
                           gene_alias       = "ALIAS",
                           gene_symbol      = "SYMBOL",
                           ensembl_gene_id  = "ENSEMBL",
                           refseq_id        = "REFSEQ",
                           "UNKNOWN"
            )
            if(prefix=="UNKNOWN"){
                stop("Currently only the following type of IDs are supported: 
entrez_id, gene_alias, ensembl_gene_id, refseq_id and gene_symbol!",
                     call.=FALSE)
            }
            tryCatch(env<-get(paste(orgAnn,prefix,"2EG",sep="")),
                     error = function(e){
                         stop(paste("Annotation database ",
                                    orgAnn,
                                    "2EG does not exist!\n
                                    \tPlease try to load annotation 
                                    database by library(",
                                    orgAnn,".db)",sep=""),call.=FALSE)
            })
            entrez <- AnnotationDbi::mget(feature_ids,env,ifnotfound=NA)
            gene_ids <- names(entrez)
            m_ent <- do.call(rbind,lapply(gene_ids,function(.ele){
                r = entrez[[.ele]]
                if(!is.na(r[1])) cbind(rep(.ele,length(r)),r)
                else {
                    if(!silence) message(paste("entrez id for '", 
                                               .ele, "' not found\n", 
                                               sep = ""))
                    c(.ele, NA)
                }
            }))
            m_ent<-as.data.frame(m_ent, stringsAsFactors=FALSE)
            m_ent<-m_ent[!is.na(m_ent[,1]), , drop=FALSE]
            colnames(m_ent)<-c(feature_id_type, "entrez_id")
            }
        entrezIDs<-as.character(m_ent$entrez_id)
        entrezIDs<-unique(entrezIDs)
        entrezIDs<-entrezIDs[!is.na(entrezIDs)]
        if(length(entrezIDs)==0){
            stop("No entrez identifier can be mapped by input data based on 
                 the feature_id_type.\nPlease consider to use correct 
                 feature_id_type, orgAnn or annotatedPeak\n",
                 call.=FALSE);
        }
        IDs2Add <- unique(IDs2Add)
        IDs2Add <- IDs2Add[IDs2Add!=feature_id_type]
        IDs <- unique(entrezIDs[!is.na(entrezIDs)])
        for(IDtoAdd in IDs2Add){
            x<-NULL
            if(!silence) message(paste("Adding",IDtoAdd,"... "))
            if(IDtoAdd!="entrez_id"){
                orgDB<-NULL
                tryCatch(orgDB<-get(paste(orgAnn,toupper(IDtoAdd),sep="")),
                         error = function(e){
                             if(!silence){
                                 message(paste("The IDs2Add you input, \"", 
                                               IDtoAdd, 
                                               "\", is not supported!\n",
                                               sep=""))
                    }
                })
                if(is.null(orgDB)){
                    IDs2Add<-IDs2Add[IDs2Add!=IDtoAdd]
                    next
                }
                if(class(orgDB)!="AnnDbBimap" & class(orgDB)!="IpiAnnDbMap"){
                    if(!silence){
                        message(paste("The IDs2Add you input, \"", 
                                      IDtoAdd, 
                                      "\", is not supported!\n",
                                      sep=""))
                    }
                    IDs2Add<-IDs2Add[IDs2Add!=IDtoAdd]
                    next
                }
                x <- AnnotationDbi::mget(IDs, orgDB,ifnotfound=NA)
                x <- sapply(x,base::paste,collapse=";")
                x <- as.data.frame(x, stringsAsFactors=FALSE)
                m_ent <- merge(m_ent, x, 
                               by.x="entrez_id",
                               by.y="row.names",
                               all.x=TRUE)
                colnames(m_ent)[length(colnames(m_ent))]<-IDtoAdd
            }
            if(!silence) message("done\n")
        }
        m_ent<-m_ent[, c(feature_id_type,IDs2Add), drop=FALSE]
}else{
    if(missing(mart) || class(mart) !="Mart"){
        stop('No valid mart object is passed in!',call.=FALSE)
    }
    IDs2Add<-unique(IDs2Add)
    IDs2Add<-IDs2Add[IDs2Add!=feature_id_type]
    tryCatch(m_ent<-
                 getBM(attributes=c(feature_id_type,IDs2Add),
                       filters = feature_id_type, 
                       values = feature_ids, mart=mart),
             error = function(e){
        stop(paste("Get error when calling getBM:", e, sep="\n"),
             call.=FALSE)
    })
    if(any(colnames(m_ent)!=c(feature_id_type, IDs2Add))) 
        colnames(m_ent) <- c(feature_id_type, IDs2Add)
}
if(!silence) message("prepare output ... ")
#dealing with multiple entrez_id for single feature_id
if(ncol(m_ent)==1) stop("None of IDs could be appended. Please double check IDs2Add.")
duplicated_ids <-
    m_ent[duplicated(m_ent[,feature_id_type]), feature_id_type]
if(length(duplicated_ids)>0){
    m_ent.duplicated <- m_ent[m_ent[,feature_id_type] %in% duplicated_ids,]
    m_ent.duplicated <- condenseMatrixByColnames(as.matrix(m_ent.duplicated),
                                                 feature_id_type)
    m_ent<-m_ent[!(m_ent[,feature_id_type] %in% duplicated_ids),]
    m_ent<-rbind(m_ent,m_ent.duplicated)
}
if (class(annotatedPeak) %in% c("GRanges")){
    #rearrange m_ent by annotatedPeak$feature
    #data.frame is very important for order...
    orderlist <- data.frame(annotatedPeak$feature)
    orderlist <- cbind(1:nrow(orderlist),orderlist)
    colnames(orderlist) <- c("orderid___",feature_id_type)
    m_ent <- merge(orderlist, m_ent, by=feature_id_type, all.x=TRUE)
    m_ent <- m_ent[order(m_ent[,"orderid___"]), 
                   c(feature_id_type, IDs2Add)]
    for(IDtoAdd in IDs2Add){
        mcols(annotatedPeak)[,IDtoAdd] <- m_ent[,IDtoAdd]
    }
}else{
    annotatedPeak <- m_ent
}
if(!silence) message("done\n")
annotatedPeak
}
