getAnnotation <-
function(mart, featureType=c("TSS","miRNA", "Exon"))
{
		featureType = match.arg(featureType)
		if (featureType == "TSS")
		{
			TSS = getBM(c("ensembl_gene_id","chromosome_name", "start_position", "end_position", "strand", "description"), mart = mart)
		}
		else if (featureType == "miRNA")
		{
			TSS = getBM(c("ensembl_gene_id","chromosome_name", "start_position", "end_position",  "strand", "description", "gene_biotype"), filter=c("biotype"),values="miRNA", mart = mart)
		}
		else if (featureType == "Exon")
		{
			TSS = getBM(c("ensembl_exon_id","chromosome_name", "exon_chrom_start", "exon_chrom_end", "strand", "description"), mart = mart)
		}
		
		TSS.unique = unique(TSS[,1:6])
		TSS.ordered <- TSS.unique[order(TSS.unique[,3]),]
		RangedData(IRanges(start=as.numeric(TSS.ordered[,3]), end=as.numeric(TSS.ordered[,4]), names= TSS.ordered[,1]),
          strand = TSS.ordered[,5],   description=TSS.ordered[,6], space = TSS.ordered[,2])
}

