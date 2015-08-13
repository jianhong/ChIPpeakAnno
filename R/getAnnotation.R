getAnnotation <-
  function(mart, featureType=c("TSS","miRNA", "Exon", "5utr", "3utr", "ExonPlusUtr", "transcript"), output=c("GRanges", "RangedData"))
  {
    featureType = match.arg(featureType)
    output = match.arg(output)
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
      TSS = getBM(c("ensembl_gene_id","chromosome_name", "start_position", "end_position",  "strand", "description", "gene_biotype"), filters=c("biotype"),values="miRNA", mart = mart)
    }
    else if (featureType == "Exon")
    {
      TSS = getBM(c("ensembl_exon_id","chromosome_name", "exon_chrom_start", "exon_chrom_end", "strand", "description"), mart = mart)
    }
    else if (featureType == "5utr")
    {
      TSS = getBM(c("ensembl_transcript_id", "chromosome_name", "5_utr_start", "5_utr_end", "strand", "description"), mart = mart)
      TSS = TSS[!is.na(TSS[,3]),]
    }
    else if (featureType == "3utr")
    {
      TSS = getBM(c("ensembl_transcript_id", "chromosome_name", "3_utr_start", "3_utr_end", "strand", "description"), mart = mart)
      
      TSS = TSS[!is.na(TSS[,3]),]
    }
    else if (featureType == "ExonPlusUtr")
    {
      TSS = getBM(c('ensembl_exon_id','chromosome_name', 'exon_chrom_start','exon_chrom_end', 
                    'strand', 'ensembl_gene_id', '5_utr_start', '5_utr_end','3_utr_start','3_utr_end', 'description'), mart=mart)
    }
    else if (featureType == "transcript")
    {
      TSS =getBM(attributes = c("ensembl_transcript_id", "chromosome_name", "transcript_start", "transcript_end", 'strand', 'description', 'ensembl_gene_id'), mart=mart)
    }
    TSS = unique(TSS)
    TSS = TSS[order(TSS[,3]),]		
    duplicatedID = TSS[duplicated(TSS[,1]),1]
    TSS = TSS[!duplicated(TSS[,1]),]
    if (length(duplicatedID) >0 && featureType != "3utr" && featureType != "5utr" && featureType != "ExonPlusUtr")
    {
      warning("Following duplicated IDs found, only one of entries of the duplicated id will be returned!")
      warning(duplicatedID)
    }
    if (featureType == "ExonPlusUtr")
    {
      re <- RangedData(IRanges(start=as.numeric(as.character(TSS[,3])),
                         end = as.numeric(as.character(TSS[,4])),
                         names= paste(as.character(TSS[,1]), as.character(TSS[,7]), as.character(TSS[,9]), sep="_")),
                 ensembl_exon_id = as.character(TSS[,1]),
                 ensembl_gene_id = as.character(TSS[,6]),
                 utr5start =  as.numeric(as.character(TSS[,7])),
                 utr5end =  as.numeric(as.character(TSS[,8])),
                 utr3start =  as.numeric(as.character(TSS[,9])),
                 utr3end =  as.numeric(as.character(TSS[,10])),
                 strand = TSS[,5],
                 description=as.character(TSS[,11]),
                 space = as.character(TSS[,2]))
    }
    else if (featureType == "transcript")
    {
      re <- RangedData(IRanges(start=as.numeric(TSS[,3]), end=as.numeric(TSS[,4]), names= as.character(TSS[,1])), strand = TSS[,5],   description=as.character(TSS[,6]), space = as.character(TSS[,2]), ensembl_gene_id=as.character(TSS[,7]))
    }
    else if (featureType == "5utr" || featureType == "3utr")
    {
      re <- RangedData(IRanges(start=as.numeric(TSS[,3]), end=as.numeric(TSS[,4]), names= paste(as.character(TSS[,1]), as.character(TSS[,3]), sep="_")), strand = TSS[,5],   description=as.character(TSS[,6]), space = as.character(TSS[,2]), ensembl_transcript_id = as.character(TSS[,1]))
    }
    else
    {
      re <- RangedData(IRanges(start=as.numeric(TSS[,3]), end=as.numeric(TSS[,4]), names= as.character(TSS[,1])), strand = TSS[,5],   description=as.character(TSS[,6]), space = as.character(TSS[,2]))
    }
    if(output=="GRanges") re <- toGRanges(re)
    re
  }
