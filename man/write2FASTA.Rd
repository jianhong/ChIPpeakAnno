\name{write2FASTA}
\alias{write2FASTA}
\title{ Write sequences to a file in fasta format}
\description{
  Write the sequences obtained from getAllPeakSequence to a file in fasta 
  format leveraging writeFASTA in Biostrings package. FASTA is a simple file 
  format for biological sequence data. A FASTA format file contains one or
  more sequences and there is a header line which begins with a > proceeding 
  each sequence.
}
\usage{
write2FASTA(mySeq, file="", width=80)
}
\arguments{
  \item{mySeq}{GRanges with varibles name and sequence ,e.g., results obtained
  from getAllPeakSequence} 
  \item{file}{Either a character string naming a file or a connection open for
  reading or writing. If "" (the default for write2FASTA), then the function 
  writes to the standard output connection (the console) unless redirected by
  sink}
  \item{width}{ The maximum number of letters per line of sequence }
}
\value{
  Output as FASTA file format to the naming file or the console. 
}
\author{ Lihua Julie Zhu}

\examples{
peaksWithSequences = GRanges(seqnames=c("1", "2"),
IRanges(start=c(1000, 2000), 
end=c(1010, 2010), 
names=c("id1", "id2")), 
sequence= c("CCCCCCCCGGGGG", "TTTTTTTAAAAAA"))

write2FASTA(peaksWithSequences, file="testseq.fasta", width=50)
}

\keyword{ misc }
