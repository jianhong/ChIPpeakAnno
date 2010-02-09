getAnnotation <-
function(mart, featureType=c("TSS","miRNA", "Exon"))
{
		featureType = match.arg(featureType)
		if (missing(mart) || class(mart) !="Mart")
		{
			stop("No valid mart object is passed in!")
		}
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
		
		TSS = unique(TSS[,1:6])
		TSS = TSS[order(TSS[,3]),]
		RangedData(IRanges(start=as.numeric(TSS[,3]), end=as.numeric(TSS[,4]), names= TSS[,1]),
          strand = TSS[,5],   description=TSS[,6], space = TSS[,2])
}

