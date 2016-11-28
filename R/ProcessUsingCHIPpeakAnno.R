#' ProcessUsingCHIPpeakAnno
#'
#' @param dir.name
#' @param input.file.pattern
#' @param out.dir.name
#'
#' @return
#' @export
#'
#' @examples
#'  dir.name="/media/H_driver/2016/Yang/MACS/MACS/"
#'  input.file.pattern="*.bed"
#'  out.dir.name="/media/H_driver/2016/Yang/Results/"
#'
#'
ProcessUsingCHIPpeakAnno <- function(dir.name,input.file.pattern,out.dir.name) {

  library(ChIPpeakAnno)
  dir.name=reformatPath(dir.name)
  out.dir.name=reformatPath(out.dir.name)

  file.name=paste0(dir.name,dir(dir.name,recursive = TRUE,pattern=input.file.pattern))
  file.name.2<-as.list(file.name)

  names(file.name.2)=sapply(strsplit(file.name,split="\\/"),"[[",8)

  print(file.name.2)

  file.name.3<-file.name.2

  sample.name<-sapply(strsplit(names(file.name.3),split="_peaks_"),"[[",1)

  names(file.name.3)=sample.name

  file.name.4 <-file.name.3[-1]

  re.out<-lapply(file.name.4,function(u){
    re=toGRanges(u,format="BED")
    #colnames(re)=c("Count","GeneName")
    re
  })

  head(re.out[[1]])

  re.out.L<-lapply(re.out,function(u){
    re=length(u)
    #colnames(re)=c("Count","GeneName")
    re
  })

  annoData <- toGRanges(EnsDb.Mmusculus.v75, feature="gene")
  annoData[1:2]

  binOverFeature(overlaps, annotationData=annoData,
                 radius=5000, nbins=20, FUN=length, errFun=0,
                 ylab="count",
                 main="Distribution of aggregated peak numbers around TSS")

  ol <- findOverlapsOfPeaks(re.out[c(2,4,1)])

  overlaps<-ol$peaklist$`11_2470IUPUI_WT_BM_SMC1_peaks.bed///13_2470IUPUI_WT_BM_Rad21_peaks.bed///10_WT_BM_ASXL1_peaks.bed`

  re<-makeVennDiagram(re.out[c(2,4,1)],NameOfPeaks=c("SMC1A", "RAD21","ASXL1"),totalTest=35000)

  #fisher exact test

  UseFisher <- function(temp.ct,index.A,index.B,totalN) {
    total.peaks=totalN
    A=sum(temp.ct[which(temp.ct[,index.A]==1&temp.ct[,index.B]==1),4])
    B=sum(temp.ct[which(temp.ct[,index.A]==1&temp.ct[,index.B]==0),4])
    C=sum(temp.ct[which(temp.ct[,index.A]==0&temp.ct[,index.B]==1),4])
    D=total.peaks-(A+B+C)
    ctb<-matrix(c(A,B,C,D),nrow = 2,dimnames =list(c("In", "Not"),c("In", "Not")))

    #re<-fisher.test(ctb)
    print(ctb)
    re.fisher<-fisher.test(ctb, alternative='greater')[c("p.value","estimate")]
    re.fisher
  }


  temp.ct<-ol$venn_cnt

  #A vs B
  index.A<-grep("SMC1",colnames(temp.ct))
  index.B<-grep("Rad21",colnames(temp.ct))
  tempRe<-UseFisher(temp.ct,index.A,index.B,35000)
  pVal.fisher.AB=tempRe$p.value
  OR.fisher.AB=tempRe$estimate

  #A vs C
  index.A<-grep("SMC1",colnames(temp.ct))
  index.B<-grep("ASXL1",colnames(temp.ct))
  tempRe<-UseFisher(temp.ct,index.A,index.B,35000)
  pVal.fisher.AC=tempRe$p.value
  OR.fisher.AC=tempRe$estimate

  #B vs C
  index.A<-grep("Rad21",colnames(temp.ct))
  index.B<-grep("ASXL1",colnames(temp.ct))
  tempRe<-UseFisher(temp.ct,index.A,index.B,35000)
  pVal.fisher.BC=tempRe$p.value
  OR.fisher.BC=tempRe$estimate

  pVal.fisher.AB
  OR.fisher.AB

  pVal.fisher.AC
  OR.fisher.AC

  pVal.fisher.BC
  OR.fisher.BC

  library(TxDb.Hsapiens.UCSC.hg19.knownGene)

  aCR<-assignChromosomeRegion(overlaps, nucleotideLevel=FALSE,
                              precedence=c("Promoters", "immediateDownstream",
                                           "fiveUTRs", "threeUTRs",
                                           "Exons", "Introns"),
                              TxDb=TxDb.Mmusculus.UCSC.mm9.knownGene)




  barplot(aCR$percentage, las=3)

  pie1(aCR$percentage,las=3)

  dc<-annoGR(TxDb.Mmusculus.UCSC.mm9.knownGene)

  #GRCm38/mm10
  dd.GRCm39.mm10<-toGRanges(EnsDb.Mmusculus.v75)
  seqinfo(dd.GRCm39.mm10)
  seqlevels(dd.GRCm39.mm10)

  seqlevels(dd.GRCm39.mm10,force=TRUE) <- c("chr1","chr10","chr11","chr12","chr13",
                        "chr14","chr15","chr16","chr17","chr18","chr19","chr2",
                        "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY")

  seqinfo(overlaps)<-seqinfo(dd.GRCm39.mm10)
  seqinfo(overlaps)

  overlaps.trimmed<-trim(overlaps, use.names=TRUE)

  library(EnsDb.Mmusculus.v79)

  #GRCm38/mm10
  dd<-toGRanges(EnsDb.Mmusculus.v79)
  seqinfo(dd)

  library(ensembldb)
  library(GenomeInfoDb)

  seqlevelsStyle(overlaps.trimmed) <- seqlevelsStyle(dd.GRCm39.mm10)

  overlaps.anno<-annoPeaks(overlaps.trimmed,dd.GRCm39.mm10)

  library(org.Mm.eg.db)

  overlaps.anno.with.entrez.id <- addGeneIDs(overlaps.anno,"org.Mm.eg.db",IDs2Add = "entrez_id")

  write.csv(as.data.frame(unname(overlaps.anno.with.entrez.id)), paste0(out.dir.name,"other_anno.csv"))

  pie1(table(overlaps.anno.with.entrez.id$insideFeature))

  library("DBI")
  over <- getEnrichedGO(overlaps.anno.with.entrez.id, orgAnn="org.Mm.eg.db",
                        maxP=.05, minGOterm=10,
                        multiAdjMethod="BH", condense=TRUE)

  # over.gene.symbol <- getEnrichedGO(overlaps.anno.with.entrez.id, orgAnn="org.Mm.eg.db",
  #                                   feature_id_type="gene_symbol",
  #                       maxP=.05, minGOterm=10,
  #                       multiAdjMethod="BH",condense=TRUE)


  head(over[["bp"]][, -c(3, 10)])

  library(org.Hs.eg.db)
  e2s = toTable(org.Mm.egSYMBOL)

  tempDS=over$bp

  tempDS2<-data.frame(apply(tempDS,1,function(u,e2c){

    #print(u[1])
    x=u[11]
    tempId<-unlist(strsplit(as.character(x),split=";"))
    index<-which(!is.na(match(e2s[,1],tempId)))
    index<-match(tempId,e2s[,1])
    geneS<-paste(e2s[index,2], collapse=";")
    geneS

    #print(geneS)

  },e2c))

  tempDS3<-cbind(tempDS,tempDS2)

  colnames(tempDS3)[12]="GeneSymbol"

  Draw4GO <- function(tempDS3) {
    x=tempDS3
    y=x[order(x$pvalue,decreasing = TRUE),]

    z=y[1:10,c(1,3,4,5,6,10)]

    Function<-z$Definition

    negative_log10p=-log10(z$pvalue)

    library(ggplot2)

    #ggplot(z, aes(x=z$go.id, y=negative_log10p,fill=factor(z$go.id)))+geom_bar(stat="identity")+geom_hline(yintercept = -log10(0.05))+coord_flip()

    ggplot(z, aes(go.id,pvalue, fill = go.id)) +
      geom_bar(stat="identity")+ scale_x_discrete(labels=z$count.InDataset, limits=factor(z$go.id))+ scale_fill_discrete(breaks = z$go.id,
      name="GO term")+theme(legend.text = element_text(colour="black", size = 11, face = "bold"))+
      theme(axis.text.x = element_text(angle = 90, hjust = 1))+ggtitle("GO enrichment analysis")+labs(x="Gene Counts",y="-log p-value")

  }


  write.table(tempDS3,file=paste0(out.dir.name,"BP_.txt"),row.names = FALSE,quote=FALSE,sep="\t")



  #anno <- annoGR(EnsDb.Hsapiens.v79)

  ree<-annotatePeakInBatch(overlaps,
                      AnnotationData=dc,
                      output="nearestBiDirectionalPromoters",
                      bindingRegion=c(-2000, 500))

  ree2 <- addGeneIDs(ree,
                              "org.Mm.eg.db",
                              IDs2Add = "entrez_id")

  re<-makeVennDiagram(re.out[c(2,4,1)],NameOfPeaks=c("SMC1A", "RAD21","ASXL1"),totalTest=35000)



  library(BSgenome.Mmusculus.UCSC.mm9)

  upseqs<-Views(Mmusculus,overlaps)

  overlaps.trimmed<-trim(overlaps, use.names=TRUE)

  mm9.S<-Seqinfo(genome="mm9")

  seqinfo(overlaps,force=TRUE) <- Seqinfo(genome="mm9")

  seqlevels(mm9.S) <- c("chr1","chr10","chr11","chr12","chr13",
                                    "chr14","chr15","chr16","chr17","chr18","chr19","chr2",
  "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY")

  seqinfo(overlaps) <- mm9.S

  goodGR <- trim(overlaps)

  overlaps.trimmed<-goodGR
  seq<-getAllPeakSequence(overlaps.trimmed,genome=Mmusculus)
  seqs.mm9 <- getSeq(Mmusculus,overlaps.trimmed)
  write2FASTA(seq, paste0(out.dir.name,"WT_triple.fa"))

  ## We can also try simulation data
  seq.sim.motif <- list(c("t", "g", "c", "a", "t", "g"),
                        c("g", "c", "a", "t", "g", "c"))
  set.seed(1)
  seq.sim <- sapply(sample(c(2, 1, 0), 1000, replace=TRUE, prob=c(0.07, 0.1, 0.83)),
                    function(x){
                      s <- sample(c("a", "c", "g", "t"),
                                  sample(100:1000, 1), replace=TRUE)
                      if(x>0){
                        si <- sample.int(length(s), 1)
                        if(si>length(s)-6) si <- length(s)-6
                        s[si:(si+5)] <- seq.sim.motif[[x]]
                      }
                      paste(s, collapse="")
                    })

  os <- oligoSummary(seq, oligoLength=6, MarkovOrder=3,
                     quickMotif=TRUE)

  zscore <- sort(os$zscore, decreasing=TRUE)
  h <- hist(zscore, breaks=100, main="Histogram of Z-score")
  text(zscore[1:2], rep(5, 2),
       labels=names(zscore[1:2]), adj=0, srt=90)

  pfms <- mapply(function(.ele, id)
    new("pfm", mat=.ele, name=paste("SAMPLE motif", id)),
    os$motifs, 1:length(os$motifs))

  motifStack(pfms[[1]])
  motifStack(pfms[[2]])
  motifStack(pfms[[3]])
  motifStack(pfms[[4]])

}
