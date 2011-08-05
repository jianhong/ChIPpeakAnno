\name{addGeneIDs}
\Rdversion{1.1}
\alias{addGeneIDs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add common IDs to annotated peaks such as gene symbol, entrez ID, ensemble gene id and refseq id.
}
\description{
Add common IDs to annotated peaks such as gene symbol, entrez ID, ensemble gene id and refseq id leveraging organism annotation dataset! For example, org.Hs.eg.db is the dataset from orgs.Hs.eg.db package for human, while org.Mm.eg.db is the dataset from the org.Mm.eg.db package for mouse}
\usage{
addGeneIDs(annotatedPeak, orgAnn, IDs2Add=c("symbol"), feature_id_type="ensembl_gene_id",silence=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{annotatedPeak}{RangedData such as data(annotatedPeak) or a vector of feature IDs}
  \item{orgAnn}{organism annotation dataset such as org.Hs.eg.db}
  \item{IDs2Add}{a vector of annotation identifiers to be added such as accnum, ensembl, ensemblprot, ensembltrans, entrez_id, enzyme, genename, pfam, pmid, prosite, refseq, symbol, unigene, and uniprot.}
  \item{feature_id_type}{type of ID: can be ensemble_gene_id, entrez_id, gene_symbol, gene_alias or refseq_id}
  \item{silence}{TRUE or FALSE. If TRUE, will not show unmapped entrez id for feature ids.}
}
\details{
	The common IDs are firstly converted to entrez IDs and then mapped to other Gene IDs.
	Parameter "IDs2Add" can be set to any combination of "accnum","ensembl","ensemblprot","ensembltrans","entrez_id","enzyme","genename","pfam","pmid","prosite","refseq","symbol","unigene","uniprot". It can be different vectors depending on the organism annotation dataset. For org.Hs.eg.db, it can also include "omim"; for org.Mm.eg.db, it can include "mgi".
	Mean of parameter "IDs2Add": 

	accnum: GenBank accession numbers;

	ensembl: Ensembl gene accession numbers;

	ensemblprot: Ensembl protein accession numbers;

	ensembltrans: Ensembl transcript accession numbers;

	entrez_id: entrez gene identifiers;

	enzyme: EC numbers;

	genename: gene name;

	pfam: Pfam identifiers;

	pmid: PubMed identifiers;

	prosite: PROSITE identifiers;

	refseq: RefSeq identifiers;

	symbol: gene abbreviations;

	unigene: UniGene cluster identifiers;

	uniprot: Uniprot accession numbers;

	omim: OMIM(Mendelian Inheritance in Man) identifiers;

	mgi: Jackson Laboratory MGI gene accession numbers
}
\value{
RangedData if the input is a RangedData or dataframe with added IDs if input is a character vector.
}
\references{
http://www.bioconductor.org/packages/release/data/annotation/
}
\author{
Jianhong Ou, Lihua Julie Zhu
}
\note{
}

\seealso{

}
\examples{
if(interactive())
{
data(annotatedPeak)
library(org.Hs.eg.db)
addGeneIDs(annotatedPeak[1:6,],"org.Hs.eg.db",c("symbol","omim"))
addGeneIDs(annotatedPeak$feature[1:6],"org.Hs.eg.db",c("symbol","genename"))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
