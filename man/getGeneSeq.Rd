\name{getGeneSeq}
\alias{getGeneSeq}
\title{ Get gene sequence using the biomaRt package  }
\description{
  Get gene sequence using the biomaRt package}
\usage{
getGeneSeq(LocationParameters, mart)
}
\arguments{
  \item{LocationParameters}{ c(ensembl_gene_id, distance from the peak to the transcription start site of the gene with the above ensemblID, upstream offset from the peak, downstream offset from the peak, Gene Start, Gene End) }
  \item{mart}{see useMart of bioMaRt package for details} 
}
\value{
  a list with the following items
  \item{feature_id}{ensemble gene ID}
  \item{distancetoFeature}{distance from the peak to the transcriptionstart site of the gene with the above ensembl gene ID}
  \item{upstream}{upstream offset from the peakStart}
  \item{downstream}{downstream offset from the peakEnd}
  \item{seq}{sequence obtained around the peak with above upstream and downstream offset}
}
\author{ Lihua Julie Zhu }
\note{ 
internal function not intended to be called directly by users
}
\examples{


if (interactive())
{
mart <- useMart(biomart="ensembl", dataset="drerio_gene_ensembl")
LocationParameters =c("ENSDARG00000054562",400, 750, 750,40454140,40454935)
getGeneSeq(LocationParameters, mart)

LocationParameters =c("ENSDARG00000054562",752, 750, 750,40454140,40454935)
getGeneSeq(LocationParameters, mart)
 
LocationParameters =c("ENSDARG00000054562",750, 750, 750,40454140,40454935)
getGeneSeq(LocationParameters, mart)
 
 LocationParameters =c("ENSDARG00000054562",-2, 750, 750,40454140,40454935)
 getGeneSeq(LocationParameters, mart)

 LocationParameters =c("ENSDARG00000054562",0, 750, 750,40454140,40454935)
 getGeneSeq(LocationParameters, mart)
 
 LocationParameters =c("ENSDARG00000054562",2, 750, 750,40454140,40454935)
 getGeneSeq(LocationParameters, mart) 
 
 LocationParameters =c("ENSDARG00000054562",1000, 750, 750,40454140,40454935)
 getGeneSeq(LocationParameters, mart)
}

}
\keyword{ internal }
\keyword{ misc }
