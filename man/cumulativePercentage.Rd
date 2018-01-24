\name{cumulativePercentage}
\alias{cumulativePercentage}
\title{Plot the cumulative percentage tag allocation in sample}
\usage{
cumulativePercentage(bamfiles, gr, input = 1, binWidth = 1000, ...)
}
\arguments{
\item{bamfiles}{Bam file names.}

\item{gr}{An object of \link[GenomicRanges:GRanges-class]{GRanges}}

\item{input}{Which file name is input. default 1.}

\item{binWidth}{The width of each bin.}

\item{...}{parameter for \link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps}.}
}
\value{
A list of data.frame with the cumulative percentages.
}
\description{
Plot the difference between the cumulative percentage tag allocation in 
paired samples.
}
\author{Jianhong Ou}
\examples{
\dontrun{
path <- system.file("extdata", "reads", package="MMDiffBamSubset")
files <- dir(path, "bam$", full.names = TRUE)
library(BSgenome.Hsapiens.UCSC.hg19)
gr <- as(seqinfo(Hsapiens)["chr1"], "GRanges")
cumulativePercentage(files, gr)
} 
 
}
\references{
Normalization, bias correction, and peak calling for ChIP-seq
Aaron Diaz, Kiyoub Park, Daniel A. Lim, Jun S. Song
Stat Appl Genet Mol Biol. Author manuscript; 
available in PMC 2012 May 3.Published in final edited form as: 
Stat Appl Genet Mol Biol. 2012 Mar 31; 11(3): 10.1515/1544-6115.1750
 /j/sagmb.2012.11.issue-3/1544-6115.1750/1544-6115.1750.xml. 
 Published online 2012 Mar 31.  doi: 10.1515/1544-6115.1750
 PMCID: PMC3342857
}
