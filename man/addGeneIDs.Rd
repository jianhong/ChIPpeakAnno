\name{addGeneIDs}
\Rdversion{1.1}
\alias{addGeneIDs}
\title{
Add common IDs to annotated peaks such as gene symbol, entrez ID, ensemble gene
id and refseq id.
}
\description{
Add common IDs to annotated peaks such as gene symbol, entrez ID, ensemble gene
id and refseq id leveraging organism annotation dataset. For example, 
org.Hs.eg.db is the dataset from orgs.Hs.eg.db package for human, while 
org.Mm.eg.db is the dataset from the org.Mm.eg.db package for mouse}
\usage{
addGeneIDs(annotatedPeak, orgAnn, IDs2Add=c("symbol"), 
           feature_id_type="ensembl_gene_id", silence=TRUE, mart)
}
\arguments{
  \item{annotatedPeak}{GRanges or a vector of feature IDs}
  \item{orgAnn}{organism annotation dataset such as org.Hs.eg.db}
  \item{IDs2Add}{a vector of annotation identifiers to be added}
  \item{feature_id_type}{type of ID to be annotated, default is 
  ensembl_gene_id}
  \item{silence}{TRUE or FALSE. If TRUE, will not show unmapped entrez id
  for feature ids.}
  \item{mart}{mart object, see \link[biomaRt:useMart]{useMart} of biomaRt 
  package for details}
}
\details{
    One of orgAnn and mart should be assigned. 
    \itemize{    
    \item If orgAnn is given, parameter feature_id_type should be 
    ensemble_gene_id, entrez_id, gene_symbol, gene_alias or refseq_id.
    And parameter IDs2Add can be set to any combination of identifiers 
    such as "accnum", "ensembl", "ensemblprot", "ensembltrans", "entrez_id",
    "enzyme", "genename", "pfam", "pmid", "prosite", "refseq", "symbol", 
    "unigene" and "uniprot". Some IDs are unique to an organism, 
    such as "omim" for org.Hs.eg.db and "mgi" for org.Mm.eg.db.

    Here is the definition of different IDs : 
        \itemize{
            \item accnum: GenBank accession numbers
            \item ensembl: Ensembl gene accession numbers
            \item ensemblprot: Ensembl protein accession numbers
            \item ensembltrans: Ensembl transcript accession numbers
            \item entrez_id: entrez gene identifiers
            \item enzyme: EC numbers
            \item genename: gene name
            \item pfam: Pfam identifiers
            \item pmid: PubMed identifiers
            \item prosite: PROSITE identifiers
            \item refseq: RefSeq identifiers
            \item symbol: gene abbreviations
            \item unigene: UniGene cluster identifiers
            \item uniprot: Uniprot accession numbers
            \item omim: OMIM(Mendelian Inheritance in Man) identifiers
            \item mgi: Jackson Laboratory MGI gene accession numbers
        }
    
    \item If mart is used instead of orgAnn, for valid parameter 
    feature_id_type and IDs2Add parameters, please refer to 
    \link[biomaRt:getBM]{getBM} in bioMart package. 
    Parameter feature_id_type should be one valid filter name listed by 
    \link[biomaRt:listFilters]{listFilters(mart)} such as ensemble_gene_id.
    And parameter IDs2Add should be one or more valid attributes name listed 
    by \link[biomaRt:listAttributes]{listAttributes(mart)} such as 
    external_gene_id, entrezgene, wikigene_name, or mirbase_transcript_name.

    }
}
\value{
GRanges if the input is a GRanges or dataframe if input is a vector.
}
\references{
http://www.bioconductor.org/packages/release/data/annotation/
}
\author{
Jianhong Ou, Lihua Julie Zhu
}

\seealso{
    \link[biomaRt:getBM]{getBM}, \link[AnnotationDbi:AnnotationDb]{AnnotationDbi}
}
\examples{
data(annotatedPeak)
library(org.Hs.eg.db)
addGeneIDs(annotatedPeak[1:6,],orgAnn="org.Hs.eg.db",
           IDs2Add=c("symbol","omim"))
addGeneIDs(annotatedPeak$feature[1:6],orgAnn="org.Hs.eg.db",
        IDs2Add=c("symbol","genename"))
if(interactive()){
mart <- useMart("ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",
                dataset="hsapiens_gene_ensembl")
##mart <- useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")
addGeneIDs(annotatedPeak[1:6,],mart=mart,IDs2Add=c("hgnc_symbol","entrezgene"))
}
}
\keyword{ misc }
