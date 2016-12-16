#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>1) {

input.file.dir=args[1]
out.file.dir=args[2]
input.file.pattern=args[3]
input.sample.file=args[4]
}

cat(input.file.dir,"\t",out.file.dir,"\t",input.file.pattern,"\n")

library(ChipSeq)

re2<-ParserBedFile4PengDiffBind(input.sample.file,input.file.dir,input.file.pattern,out.file.dir)

#BamFileSortIndexVisualization(input.file.dir,out.file.dir,genome)

#AnnotatePeak(paste0(temp3,"/"),"*macs142_peaks.bed",7,paste0(output.dir.name,"PeakAnnotation_at_",temp2),genome="Hs")

#SortBamFile(input.file.dir,input.file.pattern,as.numeric(index.file),out.file.dir,genome)

#re.from.bed.peng.4.venn<-ParserBamFile4NgsPlot(dir.name,as.numeric(index),input.file.pattern,out.dir.name,outfile)
