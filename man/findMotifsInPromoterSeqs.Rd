\name{findMotifsInPromoterSeqs}
\alias{findMotifsInPromoterSeqs}
\title{Find occurence of input motifs in the promoter regions of the input gene list}
\usage{
findMotifsInPromoterSeqs(
  patternFilePath1,
  patternFilePath2,
  findPairedMotif = FALSE,
  BSgenomeName,
  txdb,
  geneIDs,
  upstream = 5000L,
  downstream = 5000L,
  name.motif1 = "motif1",
  name.motif2 = "motf2",
  max.distance = 100L,
  min.distance = 1L,
  motif.orientation = c("both", "motif1UpstreamOfMotif2", "motif2UpstreamOfMoif1"),
  ignore.strand = FALSE,
  format = "fasta",
  skip = 0L,
  motif1LocForDistance = "end",
  motif2LocForDistance = "start",
  outfile,
  append = FALSE
)
}
\arguments{
\item{patternFilePath1}{File path containing a list of known motifs. Required}

\item{patternFilePath2}{File path containing a motif requried to be in the flanking regions of
the motif(s) in the first file, i.e, patternFilePath1. Requried if findPairedMotif is set to TRUE}

\item{BSgenomeName}{A BSgenome object. For a list of existing Bsgenomes, please refer use
the function available.genomes in BSgenome package.
For example,BSgenome.Hsapiens.UCSC.hg38 is for hg38, BSgenome.Hsapiens.UCSC.hg19 is for hg19,
BSgenome.Mmusculus.UCSC.mm10 is for mm10, BSgenome.Celegans.UCSC.ce6 is for ce6
BSgenome.Rnorvegicus.UCSC.rn5 is for rn5, BSgenome.Drerio.UCSC.danRer7 is for Zv9, and
BSgenome.Dmelanogaster.UCSC.dm3 is for dm3. Required}

\item{findPairedMotif}{Find motifs in paired configuration only or not. Default FALSE}

\item{txdb}{A TxDb object. For creating and using TxDb object, please refer
to GenomicFeatures package. For a list of existing TxDb object, please
search for annotation package starting with Txdb at
http://www.bioconductor.org/packages/release/BiocViews.html#___AnnotationData,
such as TxDb.Rnorvegicus.UCSC.rn5.refGene for rat,
TxDb.Mmusculus.UCSC.mm10.knownGene for mouse,
TxDb.Hsapiens.UCSC.hg19.knownGene and TxDb.Hsapiens.UCSC.hg38.knownGene for human,
TxDb.Dmelanogaster.UCSC.dm3.ensGene for Drosophila and
TxDb.Celegans.UCSC.ce6.ensGene for C.elegans}
\item{geneIDs}{ One or more gene entrez IDs. For example the entrez ID for EWSIR is 2130
https://www.genecards.org/cgi-bin/carddisp.pl?gene=EWSR1
You can use the addGeneIDs function in ChIPpeakAnno to convert other types of Gene IDs to entrez ID}
\item{upstream}{ Number of bases upstream of the TSS to search for the motifs. Default 5000L}
\item{downstream}{ Number of bases downstream of the TSS to search for the motifs. Default 5000L}
\item{name.motif1}{ Name of the motif in inputfilePath2 for labeling the output file column. Default motif1.
used only when searching for motifs in paired configuration}
\item{name.motif2}{ Name of the motif in inputfilePath2 for labeling the output file column. Default motif2
used only when searching for motifs in paired configuration}
\item{max.distance}{ maximum required gap between a paired motifs to be included in the output file. Default 100L}
\item{min.distance}{ Minimum required gap between a paired motifs to be included in the output file. Default 1L}
\item{motif.orientation}{ Required relative oriention between paired motifs: both means any orientation, motif1UpstreamOfMotif2 means motif1 needs to be located on the upstream of motif2, and  motif2UpstreamOfMoif1 means motif2 needs to be located on the upstream of motif1. Default both}
\item{ignore.strand}{ Specify whether paired motifs should be located on the same strand. Default FALSE}
\item{format}{ The format of the files specified in inputFilePath1 and inputFilePath2. Default fasta}
\item{skip}{ Specify number of lines to skip at the beginning of the input file. Default 0L}
\item{motif1LocForDistance}{ Specify whether to use the start or end of the motif1 location to calculate distance between paired motifs. Only applicable when findPairedMotif is set to TRUE. Default end}
\item{motif2LocForDistance}{ Specify whether to use the start or end of the motif2 location to calculate distance between paired motifs. Only applicable when findPairedMotif is set to TRUE. Default start}
\item{outfile}{ File path to save the search results}
\item{append}{ Specify whether to append the results to the specified output file, i.e., outfile.
Default FALSE}
}
\value{
A vector of numeric.
It is the background corrected log2-transformed ratios, CPMRatios or OddRatios.

An object of GRanges with metadata "tx_start", "tx_end tx_strand", "tx_id", "tx_name",
"Gene ID", and motif specific information such as motif name, motif found, motif strand etc.
}
\description{
Find occurence of input motifs in the promoter regions of the input gene list
}
\details{
This function outputs the motif occuring locations in the promoter regions of input
gene list and input motifs. It also can find paired motifs within specificed gap threshold
}
\examples{

library("BSgenome.Hsapiens.UCSC.hg38")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")

patternFilePath1 =system.file("extdata", "motifIRF4.fa", package="ChIPpeakAnno")
patternFilePath2 =system.file("extdata", "motifAP1.fa", package="ChIPpeakAnno")
pairedMotifs <- findMotifsInPromoterSeqs(patternFilePath1 = patternFilePath1,
   patternFilePath2 = patternFilePath2,
   findPairedMotif = TRUE,
   name.motif1 = "IRF4", name.motif2 = "AP1",
   BSgenomeName = BSgenome.Hsapiens.UCSC.hg38,
   geneIDs = 7486, txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
   outfile = "testPaired.xls")

unPairedMotifs <- findMotifsInPromoterSeqs(patternFilePath1 = patternFilePath1,
    BSgenomeName = BSgenome.Hsapiens.UCSC.hg38,
   geneIDs = 7486, txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
   outfile = "testUnPaired.xls")
}
\author{
Lihua Julie Zhu
}
