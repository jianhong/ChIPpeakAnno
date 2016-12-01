#! /usr/bin/env Rscript

# http://cran.r-project.org/web/packages/optparse/vignettes/optparse.pdf
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c('-p', '--PeakCall'), action='store', dest='PeakCall', type='character', default="macs2", metavar='"YYYY-MM-DD"', help='which method to call peak')
)

opt <- parse_args(OptionParser(option_list=option_list))

# print(opt$date)
cat(sprintf('select peak call method', opt$PeakCall))