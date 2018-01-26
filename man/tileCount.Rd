\name{tileCount}
\alias{tileCount}
\title{
    Perform overlap queries between reads and genome by windows
}
\description{
    tileCount extends \link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps} 
    by providing fixed window size and step to split whole genome into windows and 
    then do queries. It will return counts in each windows.
}
\usage{
tileCount(reads, genome, windowSize=1e6, step=1e6, 
                      keepPartialWindow=FALSE,
                      mode=countByOverlaps, ...)
}
\arguments{
  \item{reads}{
    A \link[GenomicRanges:GRanges-class]{GRanges}, \link[GenomicRanges:GRangesList-class]{GRangesList}
    \link[GenomicAlignments:GAlignments-class]{GAlignments}, 
    \link[GenomicAlignments:GAlignmentsList-class]{GAlignmentsList}, 
    \link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs} or 
    \link[Rsamtools:BamFile-class]{BamFileList} object that represents the data to be 
    counted by \code{\link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps}}.
  }
  \item{genome}{
    The object from/on which to get/set the sequence information.
  }
  \item{windowSize}{
    Size of windows
  }
  \item{step}{
    Step of windows
  }
  \item{keepPartialWindow}{
    Keep last partial window or not.
  }
  \item{mode}{
    mode can be one of the pre-defined count methods. see \link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps}. default is countByOverlaps, alia of countOverlaps(features, reads, ignore.strand=ignore.strand)
  }
  \item{...}{
    Additional arguments passed to \code{\link[GenomicAlignments:summarizeOverlaps-methods]{summarizeOverlaps}}. 
  }
}

\value{
A \link[SummarizedExperiment:RangedSummarizedExperiment-class]{RangedSummarizedExperiment} object. The assays slot holds the counts, rowRanges holds the annotation from genome.
}

\author{
    Jianhong Ou
}

\examples{
    fls <- list.files(system.file("extdata", package="GenomicAlignments"),
                  recursive=TRUE, pattern="*bam$", full=TRUE)
    names(fls) <- basename(fls)
    genes <- GRanges(seqlengths = c(chr2L=7000, chr2R=10000))
    se <- tileCount(fls, genes, windowSize=1000, step=500)
}
\keyword{ misc }