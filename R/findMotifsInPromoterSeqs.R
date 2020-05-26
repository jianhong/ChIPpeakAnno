#' Find occurence of input motifs in the promoter regions of the input gene list
#'
#' @details This function outputs the motif occuring locations in the promoter regions of input
#' gene list and input motifs. It also can find paired motifs within specificed gap threshold
#'
#' @param patternFilePath1 File path containing a list of known motifs. Required

#' @param patternFilePath2 File path containing a motif requried to be in the flanking regions of
#' the motif(s) in the first file, i.e, patternFilePath1. Requried if findPairedMotif is set to TRUE

#' @param findPairedMotif Find motifs in paired configuration only or not. Default FALSE

#' @param BSgenomeName A BSgenome object. For a list of existing Bsgenomes, please refer use
#' the function available.genomes in BSgenome package.
#' For example,BSgenome.Hsapiens.UCSC.hg38 is for hg38, BSgenome.Hsapiens.UCSC.hg19 is for hg19,
#' BSgenome.Mmusculus.UCSC.mm10 is for mm10, BSgenome.Celegans.UCSC.ce6 is for ce6
#' BSgenome.Rnorvegicus.UCSC.rn5 is for rn5, BSgenome.Drerio.UCSC.danRer7 is for Zv9, and
#' BSgenome.Dmelanogaster.UCSC.dm3 is for dm3. Required

#' @param txdb A TxDb object. For creating and using TxDb object, please refer
#' to GenomicFeatures package. For a list of existing TxDb object, please
#' search for annotation package starting with Txdb at
#' http://www.bioconductor.org/packages/release/BiocViews.html#___AnnotationData,
#' such as TxDb.Rnorvegicus.UCSC.rn5.refGene for rat,
#' TxDb.Mmusculus.UCSC.mm10.knownGene for mouse,
#' TxDb.Hsapiens.UCSC.hg19.knownGene and TxDb.Hsapiens.UCSC.hg38.knownGene for human,
#' TxDb.Dmelanogaster.UCSC.dm3.ensGene for Drosophila and
#' TxDb.Celegans.UCSC.ce6.ensGene for C.elegans

#' @param geneIDs One or more gene entrez IDs. For example the entrez ID for EWSIR is 2130
#' https://www.genecards.org/cgi-bin/carddisp.pl?gene=EWSR1
#' You can use the addGeneIDs function in ChIPpeakAnno to convert other types of Gene IDs to entrez ID

#' @param upstream Number of bases upstream of the TSS to search for the motifs. Default 5000L

#' @param downstream Number of bases downstream of the TSS to search for the motifs. Default 5000L

#' @param name.motif1 Name of the motif in inputfilePath2 for labeling the output file column. Default motif1.
#' used only when searching for motifs in paired configuration

#' @param name.motif2 Name of the motif in inputfilePath2 for labeling the output file column. Default motif2
#' used only when searching for motifs in paired configuration

#' @param max.distance maximum required gap between a paired motifs to be included in the output file. Default 100L

#' @param min.distance Minimum required gap between a paired motifs to be included in the output file. Default 1L

#' @param motif.orientation Required relative oriention between paired motifs: both means any orientation, motif1UpstreamOfMotif2 means motif1 needs to be located on the upstream of motif2, and  motif2UpstreamOfMoif1 means motif2 needs to be located on the upstream of motif1. Default both

#' @param ignore.strand Specify whether paired motifs should be located on the same strand. Default FALSE

#' @param format The format of the files specified in inputFilePath1 and inputFilePath2. Default fasta

#' @param skip Specify number of lines to skip at the beginning of the input file. Default 0L

#' @param motif1LocForDistance Specify whether to use the start or end of the motif1 location to calculate distance between paired motifs. Only applicable when findPairedMotif is set to TRUE. Default end

#' @param motif2LocForDistance Specify whether to use the start or end of the motif2 location to calculate distance between paired motifs. Only applicable when findPairedMotif is set to TRUE. Default start

#' @param outfile File path to save the search results

#' @param append Specify whether to append the results to the specified output file, i.e., outfile.
#' Default FALSE
#'
#' @return A vector of numeric.
#' It is the background corrected log2-transformed ratios, CPMRatios or OddRatios.
#'
#' @export
#'
#' @examples
#'
#' library("BSgenome.Hsapiens.UCSC.hg38")
#' library("TxDb.Hsapiens.UCSC.hg38.knownGene")
#'
#' patternFilePath1 =system.file("extdata", "motifIRF4.fa", package="ChIPpeakAnno")
#' patternFilePath2 =system.file("extdata", "motifAP1.fa", package="ChIPpeakAnno")
#' pairedMotifs <- findMotifsInPromoterSeqs(patternFilePath1 = patternFilePath1,
#'    patternFilePath2 = patternFilePath2,
#'    findPairedMotif = TRUE,
#'    name.motif1 = "IRF4", name.motif2 = "AP1",
#"    motif.orientation = "both",
#'    BSgenomeName = BSgenome.Hsapiens.UCSC.hg38,
#'    geneIDs = 7486, txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#'    outfile = "testPaired.xls")
#'
#' unPairedMotifs <- findMotifsInPromoterSeqs(patternFilePath1 = patternFilePath1,
#'     BSgenomeName = BSgenome.Hsapiens.UCSC.hg38,
#'    geneIDs = 7486, txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#'    outfile = "testUnPaired.xls")

#' @return An object of GRanges with metadata "tx_start", "tx_end tx_strand", "tx_id", "tx_name",
#' "Gene ID", and motif specific information such as motif name, motif found, motif strand etc.
#' @author Lihua Julie Zhu

findMotifsInPromoterSeqs <-
  function(patternFilePath1,
           patternFilePath2,
           findPairedMotif = FALSE,
           BSgenomeName,
           txdb,
           geneIDs,
           upstream = 5000L,
           downstream = 5000L,
           name.motif1 = "motif1",
           name.motif2 = "motf2",
           max.distance = 100L, min.distance = 1L,
           motif.orientation = c("both", "motif1UpstreamOfMotif2", "motif2UpstreamOfMoif1"),
           ignore.strand = FALSE,
           format = "fasta",
           skip = 0L,
           motif1LocForDistance = "end",
           motif2LocForDistance = "start",
           outfile, append = FALSE)
  {
    if (missing(patternFilePath1))
      stop("missing required parameter patternFilePath1!")
    if (missing(txdb) || !is(txdb, "TxDb"))
      stop("txdb is required as Txdb object!")
    if (missing(geneIDs))
      stop("geneIDs is required as Entrez IDs")
    motif.orientation <- match.arg(motif.orientation)
    tx <- transcriptsBy(txdb, by = "gene")
    tx <- tx[names(tx) %in% geneIDs]
    
    peaks <- promoters(tx, upstream = upstream, downstream = downstream)
    
    x1 <- do.call(rbind, lapply(1:length(peaks),function(i) {
      thisPeak <- peaks[[i]]
      mcols(thisPeak)$gene_id = names(peaks)[i]
      summarizePatternInPeaks(patternFilePath = patternFilePath1, format = format,
                              skip = skip, BSgenomeName = BSgenomeName, peaks = thisPeak)
    }
    ))
    
    
    colnames(x1)[2] = "start"
    colnames(x1)[3] = "end"
    colnames(x1)[9] <- "strand"
    colnames(x1)[11] <- "tx_start"
    colnames(x1)[12] <- "tx_end"
    colnames(x1)[14] <- "tx_strand"
    
    if (!findPairedMotif) {
      if (!missing(outfile))
        write.table(x1, file = outfile, sep ="\t", row.names = FALSE)
      toGRanges(x1)
    }
    else if(missing(patternFilePath2)) {
      stop("missing required parameter patternFilePath2!")
    }
    else {
      x2 <- do.call(rbind, lapply(1:length(peaks),function(i) {
        thisPeak <- peaks[[i]]
        mcols(thisPeak)$gene_id = names(peaks)[i]
        summarizePatternInPeaks(patternFilePath = patternFilePath2, format = format,
                                skip = skip, BSgenomeName = BSgenomeName, peaks = thisPeak)
      }
      ))
      
      colnames(x2)[2] = "start"
      colnames(x2)[3] = "end"
      colnames(x2)[9] <- "strand"
      colnames(x2)[11] <- "tx_start"
      colnames(x2)[12] <- "tx_end"
      colnames(x2)[14] <- "tx_strand"
      x1.gr <- toGRanges(x1)
      x2.gr <- toGRanges(x2)
      
      res <- annotatePeakInBatch(x1.gr, AnnotationData=x2.gr,
                                 PeakLocForDistance = motif1LocForDistance,
                                 FeatureLocForDistance = motif2LocForDistance,
                                 ignore.strand = ignore.strand)
      
      temp <- res[res$shortestDistance <= max.distance & res$shortestDistance >= min.distance,]
      
      if (motif.orientation == "motif1UpstreamMotif2" )
        temp <- temp[temp$insideFeature == "upstream",]
      if (motif.orientation == "motif1UpstreamMotif2" )
        temp <- temp[temp$insideFeature == "downstream",]
      
      temp <- as.data.frame(temp)
      
      m2 <- cbind(names(x2.gr), as.data.frame(x2.gr))
      colnames(m2)[1] <- "feature"
      m2 <- m2[, c(1:2, 7:8, 11)]
      colnames(m2)[3:5] <- paste(name.motif2, colnames(m2)[3:5])
      
      m2$seqnames <- paste("chr", m2$seqnames, sep = "" )
      
      res <- merge(m2, temp, by = c("feature", "seqnames"))
      colnames(res)[6:14] <- paste(name.motif1, colnames(res)[6:14])
      res <- res[, -c(1, 21, 28)]
      
      colnames(res)[grep("insideFeature", colnames(res))] <- paste(name.motif1, name.motif2, sep  =   "relatviePositionTo")
      
      colnames(res)[20:21] <- paste(name.motif2, colnames(res)[20:21])
      colnames(res)[22] <- paste(name.motif2, "strand", sep ="_")
      
      
      tem <- cbind(as.numeric(res[,5]), as.numeric(res[,6]), as.numeric(res[,20]), as.numeric(res[,21]))
      res2 <- cbind(seqnames = res[,1], compositMotifStart= rowMins(tem), compositMotifEnd = rowMaxs(tem), res[, -1])
      
      res2 <- res2[, -grep("Offset", colnames(res2))]
      
      d.ind <- grep("distancetoFeature", colnames(res2))
      res2 <- cbind(res2[c(1:6, d.ind, 7:dim(res2)[2])])
      res2 <- res2[, -(d.ind +1)]
      colnames(res2)[grep("distancetoFeature", colnames(res2))] <- paste(name.motif1, name.motif2, sep="distanceto")
      colnames(res2)[grep("feature_strand", colnames(res2))] <- paste(name.motif2, "strand")
      
      strand.ind <- grep("strand", colnames(res2))[1:3]
      strand <- res2[, strand.ind[1]]
      strand2 <- res2[, strand.ind[3]]
      strand <- ifelse(strand == strand2, as.character(strand), "*")
      res2 <- cbind(res2[,1:3], strand = strand, res2[,4:dim(res2)[2]])
      if (!missing(outfile)) {
        write.table(res2, file = outfile, sep ="\t", row.names = FALSE)
      }
      colnames(res2)[2:3] <- c("start", "end")
      toGRanges(res2)
    }
  }
