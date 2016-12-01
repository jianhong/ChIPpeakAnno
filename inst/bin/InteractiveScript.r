#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stats"))
# specify our desired options in a list
# by default OptionParser will add an help option equivalent to
# make_option(c("-h", "--help"), action="store_true", default=FALSE,
#               help="Show this help message and exit")

#Example: Rscript ~/ChipSeq/inst/bin/my_script.r --help

option_list <- list(
  make_option(c("-t","--orginaism"),action="store_true",default=TRUE,help="which organisim")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

cat("\n")
