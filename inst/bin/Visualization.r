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
res <- GetSampleInfo(input.sample.file,input.bam.file)

cat("please choose samples from Cell_TF :\n")

input <- file("stdin", "r")
samples.choosed <- readLines(input, n = 1)

select.sample <- as.array(samples.choosed)

print(select.sample)

ChipSeq:::configAndMultiplot(res,select.sample,output.config.dir)
