#' Obtain the TSS, exon or miRNA annotation for the specified species
#' 
#' Obtain the TSS, exon or miRNA annotation for the specified species using the
#' biomaRt package
#' 
#' 
#' @param mart A mart object, see useMart of biomaRt package for details.
#' @param featureType TSS, miRNA, Exon, 5'UTR, 3'UTR, transcript or Exon plus
#' UTR. The default is TSS.
#' @return \link[GenomicRanges:GRanges-class]{GRanges} with slot start holding
#' the start position of the feature, slot end holding the end position of the
#' feature, slot names holding the id of the feature, slot space holding the
#' chromosome location where the feature is located. In addition, the following
#' variables are included.  \item{list("strand")}{1 for positive strand and -1
#' for negative strand where the feature is located}
#' \item{list("description")}{description of the feeature such as gene}
#' @note For featureType of TSS, start is the transcription start site if
#' strand is 1 (plus strand), otherwise, end is the transcription start site.\cr\cr
#' Note that the version of the annotation db must match with the genome used for mapping because the coordinates may differ for different genome releases. For example, if you are using Mus_musculus.v103 for mapping, you'd best also use EnsDb.Mmusculus.v103 for annotation. See Examples for more info.
#' @author Lihua Julie Zhu, Jianhong Ou
#' @references Durinck S. et al. (2005) BioMart and Bioconductor: a powerful
#' link between biological biomarts and microarray data analysis.
#' Bioinformatics, 21, 3439-3440.
#' @keywords misc
#' @export
#' @import GenomicRanges
#' @importFrom biomaRt getBM
#' @examples
#' 
#' if (interactive())
#' {
#'   mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
#'   Annotation <- getAnnotation(mart, featureType="TSS")
#' }
#' 
#' Below are 3 options to fetch the annotation file.
#' Option1: with the AnnotationHub package
#' library(AnnotationHub)
#' ah <- AnnotationHub()
#' EnsDb.Mmusculus <- query(ah, pattern = c("Mus musculus", "EnsDb"))
#' EnsDb.Mmusculus
#' EnsDb.Mmusculus.v101 <- EnsDb.Mmusculus[[length(EnsDb.Mmusculus)]]
#' class(EnsDb.Mmusculus.v101)
#' 
#' Option2: with the getAnnotation() function
#' library(ChIPpeakAnno)
#' library(biomaRt)
#' listMarts()
#' mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
#' Annotation <- getAnnotation(mart)
#' Note that getAnnotation() queries biomart, which is always up-to-date.
#' 
#' Option3: build your own EnsDb package
#' This may need extra effort, and the makeEnsembldbPackage is a good starting point.


getAnnotation <- function(mart, 
                          featureType=c("TSS","miRNA", "Exon", 
                                        "5utr", "3utr", 
                                        "ExonPlusUtr", "transcript")){
    featureType = match.arg(featureType)
    if (missing(mart) || !is(mart, "Mart"))
    {
        stop("No valid mart object is passed in!")
    }
    if(featureType=="miRNA"){
        TSS = getBM(c("ensembl_gene_id","chromosome_name", "start_position", 
                      "end_position",  "strand", "description", "gene_biotype"), 
                    filters=c("biotype"),values="miRNA", mart = mart)
    }else{
        seqnames <- getBM(c("chromosome_name"), mart = mart)
        seqnames <- as.character(seqnames[, "chromosome_name"])
        TSS <- lapply(seqnames, function(.ele){
            attr <- switch(featureType, 
                           TSS=c("ensembl_gene_id", "chromosome_name", 
                                 "start_position", "end_position", 
                                 "strand", "description"),
                           Exon=c("ensembl_exon_id","chromosome_name", 
                                  "exon_chrom_start", "exon_chrom_end", 
                                  "strand", "description"),
                           "5utr"=c("ensembl_transcript_id", "chromosome_name", 
                                    "5_utr_start", "5_utr_end", "strand", 
                                    "description"),
                           "3utr"=c("ensembl_transcript_id", "chromosome_name", 
                                    "3_utr_start", "3_utr_end", "strand", 
                                    "description"),
                           ExonPlusUtr=c('ensembl_exon_id','chromosome_name', 
                                         'exon_chrom_start','exon_chrom_end', 
                                         'strand', 'ensembl_gene_id', 
                                         '5_utr_start', '5_utr_end',
                                         '3_utr_start','3_utr_end', 
                                         'description'),
                           transcript=c("ensembl_transcript_id", 
                                        "chromosome_name", 
                                        "transcript_start", "transcript_end", 
                                        'strand', 'description', 
                                        'ensembl_gene_id')
            )
            getBM(attributes=attr, filters="chromosome_name",
                  values=.ele, mart=mart)
        })
        TSS <- do.call(rbind, TSS)
    }
    
    TSS <- TSS[!is.na(TSS[, 3]), ]
    TSS = unique(TSS)
    TSS = TSS[order(TSS[,3]),]        
    duplicatedID = TSS[duplicated(TSS[,1]),1]
    TSS = TSS[!duplicated(TSS[,1]),]
    if (length(duplicatedID) >0 && 
        featureType != "3utr" && 
        featureType != "5utr" && 
        featureType != "ExonPlusUtr")
    {
        warning("Following duplicated IDs found, 
                only one of entries of the duplicated id will be returned!")
        warning(duplicatedID)
    }
    if (featureType == "ExonPlusUtr")
    {
        re <- GRanges(seqnames = as.character(TSS[,2]),
                      IRanges(start=as.numeric(as.character(TSS[,3])),
                              end = as.numeric(as.character(TSS[,4])),
                              names= make.names(as.character(TSS[,1]), 
                                                unique = TRUE)),
                      ensembl_exon_id = as.character(TSS[,1]),
                      ensembl_gene_id = as.character(TSS[,6]),
                      utr5start =  as.numeric(as.character(TSS[,7])),
                      utr5end =  as.numeric(as.character(TSS[,8])),
                      utr3start =  as.numeric(as.character(TSS[,9])),
                      utr3end =  as.numeric(as.character(TSS[,10])),
                      strand = TSS[,5],
                      description=as.character(TSS[,11]))
    }
    else if (featureType == "transcript")
    {
        re <- GRanges(seqnames = as.character(TSS[,2]),
                      IRanges(start=as.numeric(TSS[,3]), 
                              end=as.numeric(TSS[,4]), 
                              names= as.character(TSS[,1])), 
                      strand = TSS[,5],   
                      description=as.character(TSS[,6]), 
                      ensembl_gene_id=as.character(TSS[,7]))
    }
    else if (featureType == "5utr" || featureType == "3utr")
    {
        re <- GRanges(seqnames = as.character(TSS[,2]),
                      IRanges(start=as.numeric(TSS[,3]), 
                              end=as.numeric(TSS[,4]), 
                              names=make.names(as.character(TSS[,1]),
                                               unique = TRUE)), 
                      strand = TSS[,5],   
                      description=as.character(TSS[,6]), 
                      ensembl_transcript_id = as.character(TSS[,1]))
    }
    else
    {
        re <- GRanges(seqnames = as.character(TSS[,2]),
                      IRanges(start=as.numeric(TSS[,3]),
                              end=as.numeric(TSS[,4]), 
                              names= as.character(TSS[,1])), 
                      strand = TSS[,5],   
                      description=as.character(TSS[,6]))
    }
    re
}
