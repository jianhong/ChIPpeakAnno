args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>1) {
  input.sample.file=args[1]
  input.bam.file=args[2]
  output.config.dir=args[3]
}

library(ChipSeq)

ChipSeq:::configAndMultiplot(res,select.sample,output.config.dir)
