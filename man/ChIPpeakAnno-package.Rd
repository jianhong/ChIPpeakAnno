\name{ChIPpeakAnno-package}
\alias{ChIPpeakAnno-package}
\alias{ChIPpeakAnno}
\docType{package}
\title{
Batch annotation of the peaks identified from either ChIP-seq or ChIP-chip experiments.
}
\description{
 The package includes functions to retrieve the sequences around the peak, obtain enriched Gene Ontology (GO) terms, find the nearest gene, exon, miRNA or custom features such as most conserved elements and other transcription factor binding sites leveraging biomaRt, IRanges, Biostrings, BSgenome, GO.db, hypergeometric test phyper and multtest package.
}
\details{
\tabular{ll}{
Package: \tab ChIPpeakAnno\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2009-03-20\cr
License: \tab LGPL\cr
LazyLoad: \tab yes\cr
}

}
\author{
Lihua Julie Zhu, Herve Pages, Claude Gazin, Nathan Lawson, Simon Lin, David Lapointe and Michael Green


Maintainer:
 Lihua Julie Zhu <julie.zhu@umassmed.edu>
}
\references{
1. Y. Benjamini and Y. Hochberg (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. J. R. Statist. Soc. B. Vol. 57: 289-300. \cr
2. Y. Benjamini and D. Yekutieli (2001). The control of the false discovery rate in multiple hypothesis testing under dependency. Annals of Statistics. Accepted. \cr
3.  S. Durinck et al. (2005) BioMart and Bioconductor: a powerful link between biological biomarts and microarray data analysis. Bioinformatics, 21, 3439-3440. \cr
4. S. Dudoit, J. P. Shaffer, and J. C. Boldrick (Submitted). Multiple hypothesis testing in microarray experiments. \cr
5. Y. Ge, S. Dudoit, and T. P. Speed. Resampling-based multiple testing for microarray data hypothesis, Technical Report #633 of UCB Stat. http://www.stat.berkeley.edu/~gyc \cr
6. Y. Hochberg (1988). A sharper Bonferroni procedure for multiple tests of significance, Biometrika. Vol. 75: 800-802. \cr
7. S. Holm (1979). A simple sequentially rejective multiple test procedure. Scand. J. Statist.. Vol. 6: 65-70. \cr
8. N. L. Johnson,S. Kotz and A. W. Kemp (1992) Univariate Discrete Distributions, Second Edition. New York: Wiley \cr
}

\keyword{ package }
\seealso{getAnnotation, annotatePeakInBatch,  getAllPeakSequence, write2FASTA, convert2EntrezID, addAncestors, getEnrichedGO

}
\examples{

if (interactive())
{
	data(myPeakList)
        data(TSS.human.NCBI36)
	
	myPeakList1 = myPeakList[1:6,]

	annotatedPeak = annotatePeakInBatch(myPeakList1, AnnotationData=TSS.human.NCBI36)

	peaks = RangedData(IRanges(start=c(100, 500), end=c(300, 600), names=c("peak1", "peak2")), space=c("NC_008253", "NC_010468"))
library(BSgenome.Ecoli.NCBI.20080805)
	
	peaksWithSequences = getAllPeakSequence(peaks, upstream = 20,	downstream = 20, genome = Ecoli)

	write2FASTA(peaksWithSequences, file="testseq.fasta", width=50)

	library(org.Hs.eg.db)
	enrichedGO = getEnrichedGO(annotatedPeak, orgAnn ="org.Hs.eg.db", maxP=0.01, multiAdj=FALSE, minGOterm=10, multiAdjMethod="")

	enriched.biologicalprocess = enrichedGO$bp
	enriched.molecularfunction = enrichedGO$mf
	enriched.cellularcomponent = enrichedGO$cc

	data(annotatedPeak)
	y = annotatedPeak$distancetoFeature[!is.na(annotatedPeak$distancetoFeature)]
	hist(y, xlab="Distance To Nearest TSS", main="", breaks=1000, xlim=c(min(y)-100, max(y)+100))
}

}
