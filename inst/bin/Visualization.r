args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>1) {
  input.sample.file=args[1]
  input.bam.file=args[2]
  output.config.dir=args[3]
  select.sample=args[4]
}

library(ChipSeq)

res <- GetSampleInfo(input.sample.file,input.bam.file)

#select.sample <- "MDA MB 231-DD-1_cJun,MDA MB 231-1_cJun,1833-1_cJun"

tmp <- strsplit(select.sample,",")

tmp2 <- c(tmp[[1]])

ChipSeq:::configAndMultiplot(res,tmp2,output.config.dir)
