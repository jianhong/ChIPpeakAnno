#!/usr/bin/env Rscript

# source("https://bioconductor.org/biocLite.R")
# biocLite("BSgenome.Mmusculus.UCSC.mm10")

packages <- c("optparse")

null<-lapply(packages, function(xxx) suppressMessages(require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

option_list = list(make_option(c("-i", "--input_file_dir"), type="character", default=NULL,
                               help="Input file dir", metavar="character"),
                   make_option(c("-p", "--input_file_pattern"), type="character", default=NULL,
                               help="Input file dir", metavar="character"),
                   make_option(c("-o", "--back_up_dir"), type="character", default="~/Aimin/DropboxUmass/Aimin/Project",
                               help="back up dir[default= %default]", metavar="character")
);

example.use <- 'Example: Rscript ~/Aimin/DropboxUmass/Aimin/Project/Script/ChipSeqUMASS/inst/bin/Rscript/backUpFiles2DropBox.r -i ~/log -p "*html$" -o ~/Aimin/DropboxUmass/Aimin/Project'

opt_parser = OptionParser(option_list=option_list,epilogue=example.use);
opt = parse_args(opt_parser);

if (is.null(opt$input_file_dir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file dir)", call.=FALSE)
}

input.file.dir <- opt$input_file_dir
input.file.pattern <- opt$input_file_pattern
out.dir <- opt$back_up_dir

if(!dir.exists(out.dir)){dir.create(out.dir,recursive = TRUE)}

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(dirname(dirname(dirname(script.name))))
other.name <- file.path(script.basename, "R/ChipSeq.R")
print(paste("Sourcing",other.name,"from",script.name))
source(other.name)

print(input.file.dir)
print(input.file.pattern)
print(out.dir)

file.name <- list.files(input.file.dir,pattern=input.file.pattern,all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)

null <- lapply(file.name, function(u){
  x <- basename(dirname(u))
  xx <- file.path(out.dir,x)
  if(!dir.exists(xx)){dir.create(xx,recursive = TRUE)}
  cmd = paste("cd",paste0(xx,";ln -s"),paste0(u,";cd ~"),collapse = " ")
  cat(cmd,"\n")
  system(cmd)
})
