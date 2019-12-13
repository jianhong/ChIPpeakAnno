\name{getGO}
\alias{getGO}
\title{ Obtain gene ontology (GO) terms for given genes }
\description{
  Obtain gene ontology (GO) terms useing GO gene mapping package such as 
  org.Hs.db.eg to obtain the GO annotation.
}
\usage{
getGO(all.genes, orgAnn = "org.Hs.eg.db", writeTo, ID_type="gene_symbol")
}
\arguments{
  \item{all.genes}{A character vector of feature IDs}
  \item{orgAnn}{Organism annotation package such as org.Hs.eg.db for human and 
  org.Mm.eg.db for mouse, org.Dm.eg.db for fly, org.Rn.eg.db for rat, 
  org.Sc.eg.db for yeast and org.Dr.eg.db for zebrafish }
 \item{ID_type}{The feature type in annotatedPeak such as 
 ensembl_gene_id, refseq_id, gene_symbol}
  \item{writeTo}{ File path for output table}
}
\value{
An invisible table with genes and GO terms.
}

\author{ Lihua Julie Zhu }
\seealso{ getEnrichedGO }
\examples{
  if (interactive()) {
     data(annotatedPeak)
     library(org.Hs.eg.db)
     getGO(annotatedPeak[1:6]$feature, 
          orgAnn="org.Hs.eg.db", 
          ID_type="ensembl_gene_id")
}
}
\keyword{ misc }
