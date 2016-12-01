#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stats"))
# specify our desired options in a list
# by default OptionParser will add an help option equivalent to
# make_option(c("-h", "--help"), action="store_true", default=FALSE,
#               help="Show this help message and exit")

#Example: Rscript ~/ChipSeq/inst/bin/my_script.r --generator "runif" --mean 0.5 --sd 0.1

option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("-i", "--input.file.dir"), type="character", default="",
              help="Files to read in [default %default]",
              metavar="number"),
  make_option("-o", "--output.file.dir", default="",
              help = "Results to output [default \%default\"]"),
  make_option("-t","--orginaism", default="Hs",help="which organisim")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))


cat("\n")