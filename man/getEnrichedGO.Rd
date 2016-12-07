\name{getEnrichedGO}
\alias{getEnrichedGO}
\title{ Obtain enriched gene ontology (GO) terms that near the peaks }
\description{
  Obtain enriched gene ontology (GO) terms based on the features near the 
  enriched peaks using GO.db package and GO gene mapping package such as 
  org.Hs.db.eg to obtain the GO annotation and using hypergeometric test 
  (phyper) and  multtest package for adjusting p-values
}
\usage{
getEnrichedGO(annotatedPeak, orgAnn, feature_id_type="ensembl_gene_id", 
maxP=0.01, minGOterm=10, multiAdjMethod=NULL, condense=FALSE)
}
\arguments{
  \item{annotatedPeak}{A GRanges object or a vector of feature IDs}
  \item{orgAnn}{Organism annotation package such as org.Hs.eg.db for human and 
  org.Mm.eg.db for mouse, org.Dm.eg.db for fly, org.Rn.eg.db for rat, 
  org.Sc.eg.db for yeast and org.Dr.eg.db for zebrafish }
 \item{feature_id_type}{The feature type in annotatedPeak such as 
 ensembl_gene_id, refseq_id, gene_symbol or entrez_id}
  \item{maxP}{ The maximum p-value to be considered to be significant}
  \item{minGOterm}{ The minimum count in a genome for a GO term to be included}
  \item{multiAdjMethod}{ The multiple testing procedures, for details, see 
  mt.rawp2adjp in multtest package}
  \item{condense}{
    condense the results or not.
  }
}
\value{
A list with 3 elements 
 \item{\code{bp}}{ enriched biological process with the following 9 variables

    go.id:GO biological process id
    
    go.term:GO biological process term
    
    go.Definition:GO biological process description
    
    Ontology: Ontology branch, i.e. BP for biological process

    count.InDataset: count of this GO term in this dataset

    count.InGenome: count of this GO term in the genome

    pvalue: pvalue from the hypergeometric test

    totaltermInDataset: count of all GO terms in this dataset

    totaltermInGenome: count of all GO terms in the genome

}
 \item{\code{mf}}{enriched molecular function with  the following 9 variables

    go.id:GO molecular function id
    
    go.term:GO molecular function term
    
    go.Definition:GO molecular function description
    
    Ontology: Ontology branch, i.e. MF for molecular function

    count.InDataset: count of this GO term in this dataset

    count.InGenome: count of this GO term in the genome

    pvalue: pvalue from the hypergeometric test

    totaltermInDataset: count of all GO terms in this dataset

    totaltermInGenome: count of all GO terms in the genome

}
 \item{\code{cc}}{enriched cellular component  the following 9 variables

    go.id:GO cellular component id
    
    go.term:GO cellular component term
    
    go.Definition:GO cellular component description
    
    Ontology: Ontology type, i.e. CC for cellular component

    count.InDataset: count of this GO term in this dataset

    count.InGenome: count of this GO term in the genome

    pvalue: pvalue from the hypergeometric test

    totaltermInDataset: count of all GO terms in this dataset

    totaltermInGenome: count of all GO terms in the genome

}
}
\references{ Johnson, N. L., Kotz, S., and Kemp, A. W. (1992) 
Univariate Discrete Distributions, Second Edition. New York: Wiley }
\author{ Lihua Julie Zhu }
\seealso{ phyper, hyperGtest}
\examples{
  data(enrichedGO)
  enrichedGO$mf[1:10,]
  enrichedGO$bp[1:10,]
  enrichedGO$cc
  if (interactive()) {
     data(annotatedPeak)
     library(org.Hs.eg.db)
     enriched.GO = getEnrichedGO(annotatedPeak[1:6,], 
                                 orgAnn="org.Hs.eg.db", 
                                 maxP=0.01,
                                 minGOterm=10,
                                 multiAdjMethod= NULL)
     dim(enriched.GO$mf)
     colnames(enriched.GO$mf)
     dim(enriched.GO$bp)
     enriched.GO$cc
}
}
\keyword{ misc }
