\name{getEnrichedPATH}
\alias{getEnrichedPATH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Obtain enriched PATH that near the peaks }
\description{
  Obtain enriched PATH that are near the peaks using path package such as reactome.db and path mapping package such as org.Hs.db.eg to obtain the path annotation and using hypergeometric test (phyper) and  multtest package for adjusting p-values
}
\usage{
getEnrichedPATH(annotatedPeak, orgAnn, pathAnn, feature_id_type="ensembl_gene_id", 
maxP=0.01, minPATHterm=10, multiAdjMethod=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{annotatedPeak}{RangedData or GRanges such as data(annotatedPeak) or a vector of feature IDs }
  \item{orgAnn}{ organism annotation package such as org.Hs.eg.db for human and org.Mm.eg.db for mouse, org.Dm.eg.db for fly, org.Rn.eg.db for rat, org.Sc.eg.db for yeast and org.Dr.eg.db for zebrafish }
  \item{pathAnn}{ pathway annotation package such as KEGG.db, reactome.db }
  \item{feature_id_type}{the feature type in annotatedPeakRanges such as ensembl_gene_id, refseq_id, gene_symbol or entrez_id}
  \item{maxP}{ maximum p-value to be considered to be significant}
  \item{minPATHterm}{ minimum count in a genome for a path to be included  }
  \item{multiAdjMethod}{ multiple testing procedures, for details, see mt.rawp2adjp in multtest package}
}
\value{
A dataframe of enriched path with the following variables.
 \item{path.id}{KEGG PATH ID}
 \item{EntrezID}{Entrez ID}
 \item{count.InDataset}{count of this PATH in this dataset}
 \item{count.InGenome}{count of this PATH in the genome}
 \item{pvalue}{pvalue from the hypergeometric test}
 \item{totaltermInDataset}{count of all PATH in this dataset}
 \item{totaltermInGenome}{count of all PATH in the genome}
 \item{PATH}{PATH name}
}
\references{ Johnson, N. L., Kotz, S., and Kemp, A. W. (1992) Univariate Discrete Distributions, Second Edition. New York: Wiley }
\author{ Jianhong Ou }
\seealso{ phyper, hyperGtest}
\examples{
if (interactive()) {
data(annotatedPeak)
library(org.Hs.eg.db)
library(reactome.db)
enriched.PATH = getEnrichedPATH(annotatedPeak, orgAnn="org.Hs.eg.db", 
                 pathAnn="reactome.db", maxP=0.01,
                 minPATHterm=10, multiAdjMethod=NULL)
 head(enriched.PATH)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
