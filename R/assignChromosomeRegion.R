#import GenomicFeatures
## Jianhong Ou @ Mar.20, 2013
assignChromosomeRegion <-
    function(peaks.RD, exon, TSS, utr5, utr3, 
             proximal.promoter.cutoff=1000L, 
             immediate.downstream.cutoff=1000L, 
             nucleotideLevel=FALSE, 
             precedence=NULL, TxDb=NULL)
    {
        ##check inputs
        if(!is.null(TxDb)){
            if(!inherits(TxDb, "TxDb")) 
                stop("TxDb must be an object of TxDb, 
                     try\n?TxDb\tto see more info.")
            if(!inherits(peaks.RD, c("RangedData","GRanges"))) 
                stop("peaks.RD must be a GRanges object.")
            if(!is.null(precedence)) {
                if(!all(precedence %in% c("Exons", "Introns", "fiveUTRs", 
                                          "threeUTRs", "Promoters", 
                                          "immediateDownstream"))) 
                    stop("precedence must be a combination of 
                         Exons, Introns, fiveUTRs, threeUTRs, 
                         Promoters, immediateDownstream")
            }
            if(inherits(peaks.RD, "RangedData")) 
                peaks.RD <- as(peaks.RD, "GRanges")
            ignore.strand <- all(as.character(strand(peaks.RD))=="*")
            exons <- exons(TxDb, columns=NULL)
            introns <- unique(unlist(intronsByTranscript(TxDb)))
            fiveUTRs <- unique(unlist(fiveUTRsByTranscript(TxDb)))
            threeUTRs <- unique(unlist(threeUTRsByTranscript(TxDb)))
            transcripts <- unique(transcripts(TxDb, columns=NULL))
            options(warn = -1)
            try({
                promoters <- 
                    unique(promoters(TxDb, upstream=proximal.promoter.cutoff, 
                                     downstream=0))
                immediateDownstream <- 
                    unique(flank(transcripts, 
                                 width=immediate.downstream.cutoff, 
                                 start=FALSE, use.names=FALSE))
            })
            microRNAs <- tryCatch(microRNAs(TxDb), 
                                  error=function(e) return(NULL))
            tRNAs <- tryCatch(tRNAs(TxDb), error=function(e) return(NULL))
            options(warn = 0)
            annotation <- list(exons, introns, fiveUTRs, threeUTRs, 
                               promoters, immediateDownstream)
            if(!is.null(microRNAs)) 
                annotation <- c(annotation, "microRNAs"=microRNAs)
            if(!is.null(tRNAs)) 
                annotation <- c(annotation, "tRNAs"=tRNAs)
            annotation <- 
                lapply(annotation, function(.anno){mcols(.anno)<-NULL; .anno})
            names(annotation)[1:6] <- 
                c("Exons", "Introns", "fiveUTRs", "threeUTRs", 
                  "Promoters", "immediateDownstream")
            ###clear seqnames, the format should be chr+NUM
            peaks.RD <- formatSeqnames(peaks.RD)
            peaks.RD <- unique(peaks.RD)
            annotation <- lapply(annotation, formatSeqnames)
            annotation <- GRangesList(annotation)
            newAnno <- c(unlist(annotation))
            if(ignore.strand){
              newAnno.rd <- newAnno
              strand(newAnno.rd) <- "*"
              newAnno.rd <- reduce(trim(newAnno.rd))
              Intergenic.Region <- gaps(newAnno.rd, end=seqlengths(TxDb))
              Intergenic.Region <- 
                Intergenic.Region[strand(Intergenic.Region)=="*"]
            }else{
              newAnno.rd <- reduce(trim(newAnno))
              Intergenic.Region <- gaps(newAnno.rd, end=seqlengths(TxDb))
              Intergenic.Region <- 
                Intergenic.Region[strand(Intergenic.Region)!="*"]
            }
            if(!all(seqlevels(peaks.RD) %in% seqlevels(newAnno))){
                warning("peaks.RD has sequence levels not in TxDb.")
                sharedlevels <- 
                    intersect(seqlevels(newAnno), seqlevels(peaks.RD))
                peaks.RD <- keepSeqlevels(peaks.RD, sharedlevels, pruning.mode="coarse")
            }
            mcols(peaks.RD) <- NULL
            if(!is.null(precedence)){
                annotation <- 
                    annotation[unique(c(precedence,names(annotation)))]
            }
            ##    annotation$Intergenic.Region <- peaks.RD
            names(Intergenic.Region) <- NULL
            annotation$Intergenic.Region <- Intergenic.Region
            anno.names <- names(annotation)
            ol.anno <- findOverlaps(peaks.RD, annotation, ignore.strand=ignore.strand)
            if(nucleotideLevel){
              ## calculate Jaccard index
              jaccardIndex <- unlist(lapply(annotation, function(.ele){
                intersection <- intersect(.ele, peaks.RD, ignore.strand=ignore.strand)
                union <- union(.ele, peaks.RD, ignore.strand=ignore.strand)
                sum(as.numeric(width(intersection)))/sum(as.numeric(width(union)))
              }))
              jaccardIndex <- jaccardIndex[anno.names]
              names(jaccardIndex) <- anno.names
              jaccardIndex[is.na(jaccardIndex)] <- 0
              
              ## create a new annotations
              newAnno <- unlist(annotation)
              newAnno$source <- rep(names(annotation), lengths(annotation))
              newAnno.disjoin <- disjoin(newAnno, with.revmap=TRUE, ignore.strand=ignore.strand)
              if(!is.null(precedence)){
                revmap <- cbind(from=unlist(newAnno.disjoin$revmap), 
                                to=rep(seq_along(newAnno.disjoin), lengths(newAnno.disjoin$revmap)))
                revmap <- revmap[order(revmap[, "to"], revmap[, "from"]), , drop=FALSE]
                revmap <- revmap[!duplicated(revmap[, "to"]), , drop=FALSE]
                newAnno.disjoin$source <- newAnno[revmap[, "from"]]$source
              }else{
                revmap <- unlist(newAnno.disjoin$revmap)
                newAnno.disjoin <- rep(newAnno.disjoin, lengths(newAnno.disjoin$revmap))
                newAnno.disjoin$source <- newAnno[revmap]$source
              }
              ol.anno <- findOverlaps(peaks.RD, newAnno.disjoin, ignore.strand=ignore.strand)
              queryHits <- peaks.RD[queryHits(ol.anno)]
              subjectHits <- newAnno.disjoin[subjectHits(ol.anno)]
              totalLen <- sum(as.numeric(width(peaks.RD)))
              queryHits.list <- split(queryHits, subjectHits$source)
              lens <- unlist(lapply(queryHits.list, function(.ele) 
                sum(as.numeric(width(unique(.ele))))))
              percentage <- 100 * lens/totalLen
            }else{
                ##calculate Jaccard index
                ol.anno.splited <- split(queryHits(ol.anno),
                                         anno.names[subjectHits(ol.anno)])
                jaccardIndex <- unlist(lapply(anno.names, function(.name){
                    union <- length(annotation[[.name]]) + 
                        length(peaks.RD) - 
                        length(unique(subjectHits(findOverlaps(peaks.RD, 
                                                               annotation[[.name]], 
                                                               ignore.strand=ignore.strand))))
                    intersection <- length(ol.anno.splited[[.name]])
                    intersection/union
                }))
                names(jaccardIndex) <- anno.names
                ol.anno <- as.data.frame(ol.anno)
                ####keep the part only annotated in peaks.RD for peaks.RD
                ol.anno.splited <- split(ol.anno, ol.anno[,2])
                hasAnnoHits <- 
                  do.call(rbind, 
                          ol.anno.splited[names(ol.anno.splited)!=
                                            as.character(length(annotation))])
                hasAnnoHits <- unique(hasAnnoHits[,1])
                ol.anno <- 
                  ol.anno[!(ol.anno[,2]==length(annotation) & 
                              (ol.anno[,1] %in% hasAnnoHits)), ]    
                if(!is.null(precedence)){
                  ol.anno <- ol.anno[!duplicated(ol.anno[,1]), ]
                }
                ##calculate percentage
                subjectHits <-anno.names[ol.anno[,2]]
                counts <- table(subjectHits)
                percentage <- 100 * counts / length(peaks.RD)
            }
            len <- length(anno.names) - length(percentage)
            if(len>0) {
                tobeadd <- rep(0, len)
                names(tobeadd) <- anno.names[!anno.names %in% 
                                                 names(percentage)]
                percentage <- c(percentage, tobeadd)
            }
            percentage <- percentage[anno.names]
            return(list(percentage=percentage, jaccard=jaccardIndex))
        }else{
            message("Please try to use TxDb next time. Try\n
                    ?TxDb\tto see more info.")
            annotationList <- list(exon, TSS, utr5, utr3)
            names(annotationList) <- c("Exon", "TSS", "UTR5", "UTR3")
            status <- lapply(annotationList, function(.ele) {
                if(!inherits(.ele, "GRanges")){
                    stop("Annotation of exon, TSS, utr5, utr3 must 
                         be objects of GRanges.")
                }
            })
            if(!inherits(peaks.RD, "GRanges")) 
                stop("peaks.RD must be a GRanges object.") 
            ann.peaks <- annotatePeakInBatch(peaks.RD, AnnotationData = TSS)
            ann.peaks <- ann.peaks[!is.na(ann.peaks$distancetoFeature)]
            upstream <- 
                ann.peaks[ ann.peaks$insideFeature=="upstream" | 
                               (ann.peaks$distancetoFeature<0 & 
                                    ann.peaks$insideFeature == "overlapStart" & 
                                    abs(ann.peaks$distancetoFeature) >
                                    ann.peaks$shortestDistance ) | 
                               ann.peaks$insideFeature=="includeFeature" | 
                               (ann.peaks$distancetoFeature>=0 & 
                                    ann.peaks$insideFeature =="overlapStart" & 
                                    ann.peaks$distancetoFeature ==
                                    ann.peaks$shortestDistance)]
            
            proximal.promoter.n <- 
                length(upstream[upstream$distancetoFeature >= 
                                    -proximal.promoter.cutoff | 
                                    upstream$shortestDistance <= 
                                    proximal.promoter.cutoff])
            enhancer.n <- length(upstream) - proximal.promoter.n
            
            downstream <- ann.peaks[ann.peaks$insideFeature =="downstream"]
            immediateDownstream.n <- 
                length(downstream[downstream$distancetoFeature <= 
                                      immediate.downstream.cutoff,])
            enhancer.n <- enhancer.n + 
                dim(downstream[downstream$distancetoFeature > 
                                   immediate.downstream.cutoff,])
            
            inside.peaks <- 
                ann.peaks[ann.peaks$insideFeature =="inside" | 
                              ann.peaks$insideFeature ==
                              "overlapEnd" |  
                              (ann.peaks$insideFeature == "overlapStart" & 
                                   ann.peaks$distancetoFeature >=0 & 
                                   ann.peaks$distancetoFeature != 
                                   ann.peaks$shortestDistance) | 
                              (ann.peaks$insideFeature =="overlapStart" & 
                                   ann.peaks$distancetoFeature <0 & 
                                   abs(ann.peaks$distancetoFeature) ==
                                   ann.peaks$shortestDistance)]
            
            ann.utr5.peaks <- annotatePeakInBatch(inside.peaks, 
                                                  AnnotationData = utr5)
            
            proximal.promoter.n <- proximal.promoter.n + 
                length(ann.utr5.peaks[ann.utr5.peaks$insideFeature==
                                          "upstream"])
            
            utr5.n <- length(
                ann.utr5.peaks[ann.utr5.peaks$insideFeature %in% 
                                   c("includeFeature" , "inside") | 
                               (ann.utr5.peaks$insideFeature =="overlapStart" & 
                                    ann.utr5.peaks$distancetoFeature >=0 & 
                                    ann.utr5.peaks$distancetoFeature != 
                                    ann.utr5.peaks$shortestDistance)  | 
                                (ann.utr5.peaks$insideFeature =="overlapStart" & 
                                     ann.utr5.peaks$distancetoFeature <0 & 
                                     abs(ann.utr5.peaks$distancetoFeature)==
                                     ann.utr5.peaks$shortestDistance) | 
                                (ann.utr5.peaks$insideFeature =="overlapEnd" & 
                                     ann.utr5.peaks$strand=="+" & 
                                     abs(start(ann.utr5.peaks)-
                                             ann.utr5.peaks$end_position) >= 
                                     (end(ann.utr5.peaks)-
                                          ann.utr5.peaks$end_position)) | 
                                (ann.utr5.peaks$insideFeature =="overlapEnd" & 
                                     ann.utr5.peaks$strand=="-" & 
                                     abs(end(ann.utr5.peaks)-
                                             ann.utr5.peaks$start_position) >= 
                                     abs(start(ann.utr5.peaks)-
                                             ann.utr5.peaks$start_position ))])
            
            proximal.promoter.n <- 
                proximal.promoter.n +  
                length(
                    ann.utr5.peaks[
                        (ann.utr5.peaks$insideFeature =="overlapStart" & 
                             ann.utr5.peaks$distancetoFeature >=0 & 
                             ann.utr5.peaks$distancetoFeature == 
                             ann.utr5.peaks$shortestDistance)  | 
                            (ann.utr5.peaks$insideFeature =="overlapStart" & 
                                 ann.utr5.peaks$distancetoFeature <0 & 
                                 abs(ann.utr5.peaks$distancetoFeature) !=
                                 ann.utr5.peaks$shortestDistance)])
            
            downstream.utr5 <-
                ann.utr5.peaks[
                    ann.utr5.peaks$insideFeature =="downstream" |
                        (ann.utr5.peaks$insideFeature =="overlapEnd" & 
                             ann.utr5.peaks$strand=="+" & 
                             abs(start(ann.utr5.peaks)-
                                     ann.utr5.peaks$end_position) < 
                             (end(ann.utr5.peaks)-
                                  ann.utr5.peaks$end_position)) | 
                        (ann.utr5.peaks$insideFeature =="overlapEnd" & 
                             ann.utr5.peaks$strand=="-" & 
                             abs(end(ann.utr5.peaks)-
                                     ann.utr5.peaks$start_position) < 
                             abs(start(ann.utr5.peaks)-
                                     ann.utr5.peaks$start_position ))] 
            
            ann.utr3.peaks <- annotatePeakInBatch(downstream.utr5, 
                                                  AnnotationData = utr3)
            
            utr3.n <- 
                length(ann.utr3.peaks[ann.utr3.peaks$insideFeature %in% 
                                          c("includeFeature" , "overlapStart", 
                                            "overlapEnd", "inside")])
            
            rest.peaks <- ann.utr3.peaks[ann.utr3.peaks$insideFeature %in% 
                                             c("downstream", "upstream")]
            
            ann.rest.peaks <- annotatePeakInBatch(rest.peaks, 
                                                 AnnotationData = exon)
            
            intron.n <- length(ann.rest.peaks[ann.rest.peaks$insideFeature %in%
                                                  c("downstream", "upstream")])
            exon.n <- length(ann.rest.peaks) - intron.n
            
            total = length(peaks.RD)/100
            
            list( "Exons" =exon.n/total, 
                  "Introns"=intron.n/total, 
                  "fiveUTRs" = utr5.n/total, 
                  "threeUTRs" = utr3.n/total, 
                  "Promoters"= proximal.promoter.n/total, 
                  "immediate.Downstream" = immediateDownstream.n/total, 
                  "Intergenic.Region" = enhancer.n/total)
        }
    }
