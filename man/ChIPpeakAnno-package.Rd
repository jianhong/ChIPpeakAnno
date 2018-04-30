\name{ChIPpeakAnno-package}
\alias{ChIPpeakAnno-package}
\alias{ChIPpeakAnno}
\docType{package}
\title{
Batch annotation of the peaks identified from either ChIP-seq or
ChIP-chip experiments.
}
\description{
 The package includes functions to retrieve the sequences around the peak,
 obtain enriched Gene Ontology (GO) terms, find the nearest gene, exon,
 miRNA or custom features such as most conserved elements and other
 transcription factor binding sites leveraging biomaRt, IRanges, Biostrings,
 BSgenome, GO.db, hypergeometric test phyper and multtest package.
}
\details{
\tabular{ll}{
 Package: \tab ChIPpeakAnno\cr
 Type: \tab Package\cr
 Version: \tab 3.0.0\cr
 Date: \tab 2014-10-24\cr
 License: \tab LGPL\cr
 LazyLoad: \tab yes\cr
}

}
\author{
Lihua Julie Zhu, Jianhong Ou, Hervé Pagès, Claude Gazin, Nathan Lawson,
Simon Lin, David Lapointe and Michael Green


Maintainer:
Jianhong Ou <jianhong.ou@umassmed.edu>, Lihua Julie Zhu <julie.zhu@umassmed.edu>
}
\references{
1. Y. Benjamini and Y. Hochberg (1995). Controlling the false discovery rate:
a practical and powerful approach to multiple testing. J. R. Statist. Soc. B.
Vol. 57: 289-300. \cr
2. Y. Benjamini and D. Yekutieli (2001). The control of the false discovery rate
in multiple hypothesis testing under dependency. Annals of Statistics.
Accepted. \cr
3.  S. Durinck et al. (2005) BioMart and Bioconductor: a powerful link between
biological biomarts and microarray data analysis. Bioinformatics,
21, 3439-3440. \cr
4. S. Dudoit, J. P. Shaffer, and J. C. Boldrick (Submitted). Multiple
hypothesis testing in microarray experiments. \cr
5. Y. Ge, S. Dudoit, and T. P. Speed. Resampling-based multiple testing for
microarray data hypothesis, Technical Report #633 of UCB Stat.
http://www.stat.berkeley.edu/~gyc \cr
6. Y. Hochberg (1988). A sharper Bonferroni procedure for multiple tests of
significance, Biometrika. Vol. 75: 800-802. \cr
7. S. Holm (1979). A simple sequentially rejective multiple test procedure.
Scand. J. Statist.. Vol. 6: 65-70. \cr
8. N. L. Johnson,S. Kotz and A. W. Kemp (1992) Univariate Discrete
Distributions, Second Edition. New York: Wiley \cr
9. Zhu L.J. et al. (2010) ChIPpeakAnno: a Bioconductor package to annotate
ChIP-seq and ChIP-chip data. BMC Bioinformatics 2010,
11:237doi:10.1186/1471-2105-11-237.\cr
}

\keyword{ package }

\examples{
 if(interactive()){
  data(myPeakList)
  library(EnsDb.Hsapiens.v75)
  anno <- annoGR(EnsDb.Hsapiens.v75)
  annotatedPeak <-
    annotatePeakInBatch(myPeakList[1:6], AnnotationData=anno)
 }
}
