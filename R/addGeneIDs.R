addGeneIDs<-function(annotatedPeak, orgAnn, IDs2Add=c("symbol"), feature_id_type="ensembl_gene_id",silence=TRUE){
	if (missing(annotatedPeak))
	{
		stop("Missing required argument annotatedPeak!",call.=FALSE)	
	}
	if(missing(orgAnn)){
		stop('no annotation database selected',call.=FALSE)
	}
	if(!grepl(".eg.db",orgAnn,ignore.case=TRUE)){
		stop('Annotation database must be *.eg.db',call.=FALSE)
	}
	if (class(annotatedPeak) == "RangedData"){
		feature_ids <- unique(annotatedPeak$feature)
	}else if (class(annotatedPeak)  ==  "character"){
		feature_ids <- unique(annotatedPeak)
		annotatedPeak<-as.data.frame(feature_ids[!is.na(feature_ids)])
		colnames(annotatedPeak)[1] = feature_id_type
	}else{
		stop("annotatedPeak needs to be RangedData type with feature variable holding the feature id or a character vector holding the IDs of the features used to annotate the peaks!",call.=FALSE)
	}
	feature_ids<-feature_ids[!is.na(feature_ids)]
	if (length(feature_ids) == 0)
	{
		stop("There is no feature column in annotatedPeak or annotatedPeak has size 0!",call.=FALSE)
		
	}
	orgAnn<-sub("\\.db$","",orgAnn,ignore.case=TRUE)
	#get Entrez ID::entrezIDs
	if(feature_id_type=="entrez_id"){
		m_ent<-as.data.frame(feature_ids)
		m_ent<-cbind(m_ent,feature_ids)
		colnames(m_ent)<-c("entrez_id",feature_id_type)
	}else{
		prefix<-switch(feature_id_type,
					   gene_alias		= "ALIAS",
					   gene_symbol		= "SYMBOL",
					   ensembl_gene_id	= "ENSEMBL",
					   refseq_id		= "REFSEQ",
					   "UNKNOWN"
					   )
		if(prefix=="UNKNOWN"){
			stop("Currently only the following type of IDs are supported: gene_alias, 
				 ensembl_gene_id, refseq_id and gene_symbol!",call.=FALSE)
		}
		tryCatch(env<-get(paste(orgAnn,prefix,"2EG",sep="")),error = function(e){
					stop(paste("Annotation database ",orgAnn,"2EG does not exist!\n\tPlease try to load annotation database by library(",orgAnn,".db)",sep=""),call.=FALSE)
				 })
		#entrez <- AnnotationDbi::mget(feature_ids,env,ifnotfound=NA)
		entrez <- mget(feature_ids,env,ifnotfound=NA)
		gene_ids <- names(entrez)
		m_ent <- do.call(rbind,lapply(gene_ids,function(.ele){
									  r = entrez[names(entrez)==.ele]
									  if(!is.na(r[1])) cbind(rep(.ele,length(r[[1]])),r[[1]])
									  else {
										if(!silence) cat(paste("entrez id for '", .ele, "' not found\n", sep = ""))
									  }
									  }))
		if(is.null(m_ent)){
			stop("No entrez identifier can be mapped by input data based on the feature_id_type.\nPlease consider to use correct feature_id_type, orgAnn or annotatedPeak\n",call.=FALSE);
		}
		m_ent<-as.data.frame(m_ent)
		m_ent<-m_ent[!is.na(m_ent[,1]),]
		colnames(m_ent)<-c(feature_id_type,"entrez_id")
	}
	entrezIDs<-as.character(m_ent$entrez_id)
	IDs2Add<-unique(IDs2Add)
	IDs<- unique(entrezIDs[!is.na(entrezIDs)])
	for(IDtoAdd in IDs2Add){
		gc()
		x<-NULL
		cat("Adding",IDtoAdd,"... ")
		if(IDtoAdd=="entrez_id"){
			x<-m_ent
		}else{
			orgDB<-NULL
			tryCatch(orgDB<-get(paste(orgAnn,toupper(IDtoAdd),sep="")),error = function(e){
					 cat(paste("The IDs2Add you input, \"", IDtoAdd, "\", is not supported!\n",sep=""))
					 })
			if(is.null(orgDB)){
				next
			}
			if(class(orgDB)!="AnnDbBimap" & class(orgDB)!="IpiAnnDbMap"){
				cat(paste("The IDs2Add you input, \"", IDtoAdd, "\", is not supported!\n",sep=""))
				next
			}
			#x <- AnnotationDbi::mget(IDs, orgDB,ifnotfound=NA)
			x <- mget(IDs, orgDB,ifnotfound=NA)
			x <- lapply(x,function(.ele){
						tmp<-paste(as.character(unlist(.ele)),sep="",collapse=";")
						gsub(";NA","",tmp)
						})
			x <- do.call(rbind,x)
			x <- unique(x)
			colnames(x)<-c(IDtoAdd)
			x <- merge(m_ent,x,by.x="entrez_id",by.y="row.names")
			x <- x[,c(feature_id_type,IDtoAdd)]
		}
		#merge duplicate feature_id
		x <- as.matrix(x)
		x <- split(x[,IDtoAdd],x[,feature_id_type])
		x <- lapply(x,paste,sep="",collapse=";")
		x <- do.call(rbind,x)
		x <- cbind(rownames(x),x)
		colnames(x)<-c(feature_id_type,IDtoAdd)
		if (class(annotatedPeak) == "RangedData"){
			tobeAdd<-do.call(rbind,lapply(annotatedPeak$feature,function(.ele){
										  r=as.character(x[x[,feature_id_type]==.ele,IDtoAdd])
										  if(length(r)>0){
										  c(.ele,r[[1]][1])
										  }else{
										  c(.ele,NA)
										  }}))
			colnames(tobeAdd)<-c(feature_id_type,IDtoAdd)
			annotatedPeak$tobeAdd<-tobeAdd[,IDtoAdd]
			cname<-colnames(annotatedPeak)
			cname[length(cname)]<-IDtoAdd
			colnames(annotatedPeak)<-cname
		}else{
			annotatedPeak<-merge(annotatedPeak,x,by=feature_id_type,all.x=TRUE)
		}
		cat("done\n")
	}
	annotatedPeak
}