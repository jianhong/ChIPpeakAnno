#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>1) {

input.file.dir=args[1]
output.file.dir=args[2]
genome=args[3]
}

cat(input.file.dir,"\t",output.file.dir,"\t",genome,"\n")

library(ChipSeq)

BamFileSortIndexVisualization(input.file.dir,output.file.dir,genome)