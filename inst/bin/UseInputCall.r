#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>1) {

input.sample.file=args[1]
input.bam.file=args[2]
genome=args[3]
output.dir=args[4]
peakcaller=args[5]
peakPvalue=args[6]
}

cat(input.sample.file,"\t",input.bam.file,"\t",genome,"\t",output.dir,"\t",peakcaller,"\t",peakPvalue,"\n")

library(ChipSeq)

#PeakCallAndAnnotation(input.file.dir,out.file.dir,genome)
re<-peakcallwithinput(input.sample.file,input.bam.file,genome,output.dir,peakcaller,peakPvalue)

#BamFileSortIndexVisualization(input.file.dir,out.file.dir,genome)

#AnnotatePeak(paste0(temp3,"/"),"*macs142_peaks.bed",7,paste0(output.dir.name,"PeakAnnotation_at_",temp2),genome="Hs")

#SortBamFile(input.file.dir,input.file.pattern,as.numeric(index.file),out.file.dir,genome)

#re.from.bed.peng.4.venn<-ParserBamFile4NgsPlot(dir.name,as.numeric(index),input.file.pattern,out.dir.name,outfile)
