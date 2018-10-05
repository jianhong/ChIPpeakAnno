#!/usr/bin/env Rscript

# source("https://bioconductor.org/biocLite.R")
# biocLite("BSgenome.Mmusculus.UCSC.mm10")

packages <- c("optparse","ChIPpeakAnno","ChIPseeker","EnsDb.Mmusculus.v75","DBI","tidyverse","annotate",
              "reactome.db","dplyr","ggplot2","ggpubr","stringr","grid","Vennerable","devtools","reshape2",
              "binaryLogic","org.Hs.eg.db","org.Mm.eg.db","eulerr","DiffBind","LOLA",
              "TxDb.Hsapiens.UCSC.hg38.knownGene","EnsDb.Hsapiens.v75","BSgenome.Hsapiens.UCSC.hg19",
              "BSgenome.Mmusculus.UCSC.mm10")

null<-lapply(packages, function(xxx) suppressMessages(require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

option_list = list(make_option(c("-i", "--input_file_dir"), type="character", default=NULL,
                               help="Input file dir", metavar="character"),
                   make_option(c("-o", "--out_dir"), type="character", default=".",
                               help="output fir[default= %default]", metavar="character")
);

example.use <- 'Example: Rscript /home/ay64w/Project/ChipSeqUMASS/inst/bin/Rscript/ChipSeqUMASS_getDB.r -i /project/umw_paul_kaufman/Aimin/Results_10_6_2018_ctcf -o /project/umw_paul_kaufman/Aimin/Results_10_6_2018_ctcf'

opt_parser = OptionParser(option_list=option_list,epilogue=example.use);
opt = parse_args(opt_parser);

if (is.null(opt$input_file_dir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file dir)", call.=FALSE)
}

input.file.dir <- opt$input_file_dir
out.dir <- opt$out_dir

if(!dir.exists(out.dir)){dir.create(out.dir,recursive = TRUE)}

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(dirname(dirname(dirname(script.name))))
other.name <- file.path(script.basename, "R/ChipSeq.R")
print(paste("Sourcing",other.name,"from",script.name))
source(other.name)

print(input.file.dir)
print(out.dir)

#AnnotatePeakUMASS(input.file.dir,input.file.pattern,out.dir,genome=input.genome)
#/project/umw_paul_kaufman/Aimin/Results_10_6_2018_ctcf/
XL.nonXL.subset <- readRDS(file.path(input.file.dir,"XL_nonXL.rds"))
RegionDB.ctcf <- readRDS(file.path(input.file.dir,"ctcfRegionDB.rds"))

res <- getOverlapWithOther(XL.nonXL.subset,select.query.peak.index= NULL,RegionDB.ctcf,out.dir)


