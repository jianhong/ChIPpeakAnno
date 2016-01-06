\name{oligoSummary}
\alias{oligoSummary}
\title{
    Output a summary of consensus in the peaks
}
\description{
   Calculate the z-scores of all combinations of oligonucleotide in a given length
   by Markove chain.
}
\usage{
oligoSummary(sequence, oligoLength = 6L, freqs = NULL, 
             MarkovOrder = 3L, quickMotif = FALSE, revcomp=FALSE)
}
\arguments{
  \item{sequence}{
    The sequences packaged in DNAStringSet, DNAString object or output of 
    function \link{getAllPeakSequence}.
}
  \item{oligoLength}{
    The length of oligonucleotide.
}
  \item{freqs}{
    Output of function \link{frequency}.
}
  \item{MarkovOrder}{
    The order of Markov chain.
}
  \item{quickMotif}{
    Generate the motif by z-score of not.
}
  \item{revcomp}{
    Consider both the given strand and the reverse complement strand when searching for motifs in a complementable alphabet (ie DNA). Default, FALSE.
}
}
\value{
A list is returned.
\item{zscore}{A numeric vector. The z-scores of each oligonucleotide.}
\item{counts}{A numeric vector. The counts number of each oligonucleotide.}
\item{motifs}{a list of motif matrix.}
}
\references{
van Helden, Jacques, Marcel li del Olmo, and Jose E. Perez-Ortin. "Statistical analysis of yeast genomic downstream sequences reveals putative polyadenylation signals." Nucleic Acids Research 28.4 (2000): 1000-1010.
}
\author{
Jianhong Ou
}
\seealso{
See Also as \code{\link{frequency}}
}
\examples{
    if(interactive()){
        data(annotatedPeak)
        library(BSgenome.Hsapiens.UCSC.hg19)
        seq <- getAllPeakSequence(annotatedPeak[1:100], 
                     upstream=20, 
                     downstream=20, 
                     genome=Hsapiens)
        oligoSummary(seq)
    }
}
\keyword{ misc }
