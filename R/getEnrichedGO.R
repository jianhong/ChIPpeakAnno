getEnrichedGO <-
function(annotatedPeak, orgAnn, feature_id_type="ensembl_gene_id", maxP=0.01, multiAdj=FALSE, minGOterm=10, multiAdjMethod="")
{	
	if (missing(annotatedPeak))
	{
		stop("Missing required argument annotatedPeak!")	
	}
	if (missing(orgAnn))
	{
		message("No valid organism specific GO gene mapping package as argument orgAnn is passed in!")
		stop("Please refer http://www.bioconductor.org/packages/release/data/annotation/ for available org.xx.eg.db packages")
	}
	GOgenome = sub(".db","",orgAnn)
	if (nchar(GOgenome) <1)
	{
		message("No valid organism specific GO gene mapping package as parameter orgAnn is passed in!")
		stop("Please refer http://www.bioconductor.org/packages/release/data/annotation/ for available org.xx.eg.db packages")
	}
	if (class(annotatedPeak) == "RangedData")
	{
		feature_ids = unique(annotatedPeak$feature)
	}
	else if (class(annotatedPeak)  ==  "character")
	{
		feature_ids = unique(annotatedPeak)
	}
	else
	{
		stop("annotatedPeak needs to be RangedData type with feature variable holding the feature id or a character vector holding the IDs of the features used to annotate the peaks!")
	}
	if (feature_id_type == "entrez_id")
	{
		entrezIDs <- feature_ids
	}
	else
	{
		entrezIDs <- convert2EntrezID(feature_ids, orgAnn, feature_id_type)
	}
		
	goAnn <- get(paste(GOgenome,"GO", sep=""))
	mapped_genes <- mappedkeys(goAnn)
	xx <- as.list(goAnn[mapped_genes])
	all.GO= matrix(unlist(unlist(xx)),ncol=3,byrow=TRUE)

	this.GO <- do.call(rbind, lapply(entrezIDs,function(x1)
	{
		temp = unlist(xx[names(xx) ==x1])
		if (length(temp) >0)
		{
			matrix(temp,ncol=3,byrow=TRUE)
		}
	}))

	bp.go.this  = addAncestors(this.GO[this.GO[,3]=="BP",1],"bp")
	cc.go.this = addAncestors(this.GO[this.GO[,3]=="CC",1], "cc")
	mf.go.this = addAncestors(this.GO[this.GO[,3]=="MF",1],"mf")

	bp.go.all  = addAncestors(all.GO[all.GO[,3]=="BP",1], "bp")
	cc.go.all = addAncestors(all.GO[all.GO[,3]=="CC",1], "cc")
	mf.go.all = addAncestors(all.GO[all.GO[,3]=="MF",1], "mf")
		
	total.mf = length(mf.go.all)
	total.cc = length(cc.go.all)
	total.bp = length(bp.go.all)
	this.mf = length(mf.go.this)
	this.cc = length(cc.go.this)
	this.bp = length(bp.go.this)
	
	this.bp.count = getUniqueGOidCount(as.character(bp.go.this[bp.go.this!=""]))
	this.mf.count = getUniqueGOidCount(as.character(mf.go.this[mf.go.this!=""]))
	this.cc.count = getUniqueGOidCount(as.character(cc.go.this[cc.go.this!=""]))
	
	all.bp.count = getUniqueGOidCount(as.character(bp.go.all[bp.go.all!=""]))
	all.mf.count = getUniqueGOidCount(as.character(mf.go.all[mf.go.all!=""]))
	all.cc.count = getUniqueGOidCount(as.character(cc.go.all[cc.go.all!=""]))
	
	bp.selected = hyperGtest(all.bp.count,this.bp.count, total.bp, this.bp)
	mf.selected = hyperGtest(all.mf.count,this.mf.count, total.mf, this.mf)
	cc.selected = hyperGtest(all.cc.count,this.cc.count, total.cc, this.cc)
	
	bp.selected = data.frame(bp.selected)
	mf.selected = data.frame(mf.selected)
	cc.selected = data.frame(cc.selected)
	
	colnames(bp.selected) = c("go.id", "count.InDataset", "count.InGenome", "pvalue", "totaltermInDataset", "totaltermInGenome")
	colnames(mf.selected) = c("go.id", "count.InDataset", "count.InGenome", "pvalue", "totaltermInDataset", "totaltermInGenome")
	colnames(cc.selected) = c("go.id", "count.InDataset", "count.InGenome", "pvalue","totaltermInDataset", "totaltermInGenome")
	
	if (multiAdj == FALSE | multiAdjMethod == "")
	{
		bp.s = bp.selected[as.numeric(as.character(bp.selected[,4]))<maxP & as.numeric(as.character(bp.selected[,3]))>=minGOterm,]
		mf.s = mf.selected[as.numeric(as.character(mf.selected[,4]))<maxP & as.numeric(as.character(mf.selected[,3]))>=minGOterm,]
		cc.s = cc.selected[as.numeric(as.character(cc.selected[,4]))<maxP & as.numeric(as.character(cc.selected[,3]))>=minGOterm,]
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
		
		bp.s = bp1[as.numeric(as.character(bp1[,dim(bp1)[2]]))<maxP &  !is.na(bp1[,dim(bp1)[2]]) & as.numeric(as.character(bp1[,4]))>=minGOterm,]
		mf.s = mf1[as.numeric(as.character(mf1[,dim(mf1)[2]]))<maxP & !is.na(mf1[,dim(mf1)[2]]) & as.numeric(as.character(mf1[,4]))>=minGOterm,]
		cc.s = cc1[as.numeric(as.character(cc1[,dim(cc1)[2]]))<maxP & !is.na(cc1[,dim(cc1)[2]]) & as.numeric(as.character(cc1[,4]))>=minGOterm,]
	}
	
    xx <- as.list(GOTERM)
	goterm.bp = do.call(rbind, lapply(as.character(bp.s$go.id),function(x1)
		{
			r= xx[names(xx)==x1]
			if (length(r) >0)
			{
				c(GOID(r[[1]]),Term(r[[1]]),Definition(r[[1]]), Ontology(r[[1]]))
			}
		}))
	if (length(goterm.bp) <1)
	{
		goterm.bp =matrix(ncol=4)
	}
	colnames(goterm.bp) = c("go.id", "go.term", "Definition", "Ontology")
	
	goterm.mf = do.call(rbind, lapply(as.character(mf.s$go.id),function(x1)
		{
			r= xx[names(xx)==x1]
			if (length(r) >0)
			{
				c(GOID(r[[1]]),Term(r[[1]]),Definition(r[[1]]), Ontology(r[[1]]))
			}
		}))
	if (length(goterm.mf) <1)
	{
		goterm.mf =matrix(ncol=4)
	}
	colnames(goterm.mf) = c("go.id", "go.term", "Definition", "Ontology")
	
	goterm.cc = do.call(rbind, lapply(as.character(cc.s$go.id),function(x1)
		{
			r= xx[names(xx)==x1]
			if (length(r) >0)
			{
				c(GOID(r[[1]]),Term(r[[1]]),Definition(r[[1]]), Ontology(r[[1]]))
			}
		}))
	if (length(goterm.cc) <1)
	{
		goterm.cc =matrix(ncol=4)
	}
	colnames(goterm.cc) = c("go.id", "go.term", "Definition", "Ontology")
	
	bp.selected = merge(goterm.bp, bp.s,by="go.id")
	mf.selected = merge(goterm.mf, mf.s,by="go.id")
	cc.selected = merge(goterm.cc,cc.s,by="go.id")
	list(bp=bp.selected, mf=mf.selected, cc=cc.selected)
}
