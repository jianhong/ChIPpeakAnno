#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>1) {

input.file.dir=args[1]
#index.file=args[2]
#input.file.pattern=args[3]
out.file.dir=args[2]
genome=args[3]
}

cat(input.file.dir,"\t",out.file.dir,"\t",genome,"\n")

library(ChipSeq)

PeakCallAndAnnotation(input.file.dir,out.file.dir,genome)

#AnnotatePeak(paste0(temp3,"/"),"*macs142_peaks.bed",7,paste0(output.dir.name,"PeakAnnotation_at_",temp2),genome="Hs")

#SortBamFile(input.file.dir,input.file.pattern,as.numeric(index.file),out.file.dir,genome)

#re.from.bed.peng.4.venn<-ParserBamFile4NgsPlot(dir.name,as.numeric(index),input.file.pattern,out.dir.name,outfile)
