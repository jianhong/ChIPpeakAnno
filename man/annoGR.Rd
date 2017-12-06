\name{annoGR-class}
\docType{class}
\alias{annoGR-class}
\alias{annoGR}
\alias{annoGR,GRanges-method}
\alias{annoGR,TxDb-method}
\alias{annoGR,EnsDb-method}
\alias{info}
\alias{info,annoGR-method}
\alias{coerce,GRanges,annoGR-method}
\alias{coerce,annoGR,GRanges-method}
\title{
    Class \code{annoGR}
}
\description{
    An object of class \code{annoGR} represents the annotation data could be
    used by annotationPeakInBatch.
}
\section{Objects from the Class}{
    Objects can be created by calls of the form
    \code{new("annoGR", date, elementMetadata, feature, mdata, ranges,
                        seqinfo, seqnames, source, strand)}
}
\section{Slots}{
  \describe{
    \item{seqnames, ranges, strand, elementMetadata, seqinfo}{
    slots inherit from \link[GenomicRanges]{GRanges}. 
    The ranges must have unique names.}
    
    \item{\code{source}}{character, where the annotation comes from}
    
    \item{\code{date}}{a \link{Date} object}
    
    \item{\code{feature}}{annotation type, could be "gene", "exon", 
    "transcript", "CDS", "fiveUTR", "threeUTR", "microRNA", "tRNAs", 
    "geneModel" for \link[GenomicFeatures]{TxDb} object, or "gene", 
    "exon" "transcript" for \link[ensembldb]{EnsDb} object}
    
    \item{\code{mdata}}{data frame, metadata from annotation}
  }
}
\section{Coercion}{
 \describe{
   \item{}{
    as(from, "annoGR"): Creates a annoGR object from a GRanges object.
   }
   \item{}{
    as(from, "GRanges"): Create a GRanges object from a annoGR object.
   }
 }
}
\section{Methods}{
  \describe{
    \item{info}{Print basic info for annoGR object}
    \item{annoGR("TxDb"), annoGR("EnsDb")}{Create a annoGR object from 
    \link[GenomicFeatures]{TxDb} or \link[ensembldb]{EnsDb} object}
  }
}
\usage{
    \S4method{annoGR}{GRanges}(ranges, feature="group", date, ...)
    \S4method{annoGR}{TxDb}(ranges, feature=c(
                            "gene", "transcript", "exon",
                            "CDS", "fiveUTR", "threeUTR",
                            "microRNA", "tRNAs", "geneModel"),
                   date, source, mdata, OrganismDb)
    \S4method{annoGR}{EnsDb}(ranges, 
                   feature=c("gene", "transcript", "exon", "disjointExons"),
                   date, source, mdata)
}
\arguments{
    \item{ranges}{an object of \link[GenomicRanges]{GRanges}, 
        \link[GenomicFeatures]{TxDb} or \link[ensembldb]{EnsDb}}
    \item{feature}{annotation type}
    \item{date}{a \link{Date} object}
    \item{...}{could be following parameters}
    \item{source}{character, where the annotation comes from}
    \item{mdata}{data frame, metadata from annotation}
    \item{OrganismDb}{an object of OrganismDb. It is used 
        for extracting gene symbol for geneModel group for 
        \link[GenomicFeatures]{TxDb}}
}
\author{
    Jianhong Ou
}
\examples{
    if(Sys.getenv("USER")=="jianhongou"){
        library(EnsDb.Hsapiens.v79)
        anno <- annoGR(EnsDb.Hsapiens.v79)
    }
}
\keyword{ classes }