#' Obtain the distance to the nearest TSS, miRNA, and/or exon for a list of
#' peaks
#' 
#' Obtain the distance to the nearest TSS, miRNA, exon et al for a list of peak
#' locations leveraging IRanges and biomaRt package
#' 
#' 
#' @param myPeakList A \link[GenomicRanges:GRanges-class]{GRanges} object
#' @param mart A mart object, used if AnnotationData is not supplied, see
#' useMart of bioMaRt package for details
#' @param featureType A charcter vector used with mart argument if
#' AnnotationData is not supplied; it's value is "TSS"", "miRNA"" or "Exon"
#' @param AnnotationData A \link[GenomicRanges:GRanges-class]{GRanges} or
#' \link{annoGR} oject.  It can be obtained from function getAnnotation or
#' customized annotation of class GRanges containing additional variable:
#' strand (1 or + for plus strand and -1 or - for minus strand). Pre-compliled
#' annotations, such as TSS.human.NCBI36, TSS.mouse.NCBIM37, TSS.rat.RGSC3.4
#' and TSS.zebrafish.Zv8, are provided by this package (attach them with data()
#' function). Another method to provide annotation data is to obtain through
#' biomaRt real time by using the parameters of mart and featureType
#' @param output \describe{ \item{nearestLocation (default)}{will output the
#' nearest features calculated as PeakLoc - FeatureLocForDistance}
#' \item{overlapping}{will output overlapping features with maximum gap
#' specified as maxgap between peak range and feature range}
#' \item{shortestDistance}{will output nearest features} \item{both}{will
#' output all the nearest features, in addition, will output any features that
#' overlap the peak that is not the nearest features}
#' \item{upstream&inside}{will output all upstream and overlapping features
#' with maximum gap} \item{inside&downstream}{will output all downstream and
#' overlapping features with maximum gap} \item{upstream}{will output all
#' upstream features with maximum gap.} \item{downstream}{will output all
#' downstream features with maximum gap.} \item{upstreamORdownstream}{will
#' output all upstream features with maximum gap or downstream with maximum
#' gap} \item{nearestBiDirectionalPromoters}{will use \link{annoPeaks} to
#' annotate peaks. Nearest promoters from both direction of the peaks (strand
#' is considered). It will report bidirectional promoters if there are
#' promoters in both directions in the given region (defined by bindingRegion).
#' Otherwise, it will report the closest promoter in one direction.} }
#' @param multiple Not applicable when output is nearest. TRUE: output multiple
#' overlapping features for each peak. FALSE: output at most one overlapping
#' feature for each peak. This parameter is kept for backward compatibility,
#' please use select.
#' @param maxgap The maximum \emph{gap} that is allowed between 2 ranges for
#' the ranges to be considered as overlapping. The \emph{gap} between 2 ranges
#' is the number of positions that separate them. The \emph{gap} between 2
#' adjacent ranges is 0. By convention when one range has its start or end
#' strictly inside the other (i.e. non-disjoint ranges), the \emph{gap} is
#' considered to be -1.
#' @param PeakLocForDistance Specify the location of peak for calculating
#' distance,i.e., middle means using middle of the peak to calculate distance
#' to feature, start means using start of the peak to calculate the distance to
#' feature. To be compatible with previous version, by default using start
#' @param FeatureLocForDistance Specify the location of feature for calculating
#' distance,i.e., middle means using middle of the feature to calculate
#' distance of peak to feature, start means using start of the feature to
#' calculate the distance to feature, TSS means using start of feature when
#' feature is on plus strand and using end of feature when feature is on minus
#' strand, geneEnd means using end of feature when feature is on plus strand
#' and using start of feature when feature is on minus strand. To be compatible
#' with previous version, by default using TSS
#' @param select "all" may return multiple overlapping peaks, "first" will
#' return the first overlapping peak, "last" will return the last overlapping
#' peak and "arbitrary" will return one of the overlapping peaks.
#' @param ignore.strand When set to TRUE, the strand information is ignored in
#' the annotation.
#' @param bindingRegion Annotation range used for \link{annoPeaks}, which is a
#' vector with two integer values, default to c (-5000, 5000). The first one
#' must be no bigger than 0. And the sec ond one must be no less than 1. Once
#' bindingRegion is defined, annotation will based on \link{annoPeaks}. Here is
#' how to use it together with the parameter output and FeatureLocForDistance.
#' \itemize{ \item To obtain peaks with nearest bi-directional promoters within
#' 5kb upstream and 3kb downstream of TSS, set output =
#' "nearestBiDirectionalPromoters" and bindingRegion = c(-5000, 3000) \item To
#' obtain peaks within 5kb upstream and up to 3kb downstream of TSS within the
#' gene body, set output="overlapping", FeatureLocForDistance="TSS" and
#' bindingRegion = c(-5000, 3000) \item To obtain peaks up to 5kb upstream
#' within the gene body and 3kb downstream of gene/Exon End, set
#' output="overlapping", FeatureLocForDistance="geneEnd" and bindingRegion =
#' c(-5000, 3000) \item To obtain peaks from 5kb upstream to 3kb downstream of
#' genes/Exons, set output="overlapping", bindingType = "fullRange" and
#' bindingRegion = c(-5000, 3000) } For details, see \link{annoPeaks}.
#' @param ... Parameters could be passed to \link{annoPeaks}
#' @return An object of \link[GenomicRanges:GRanges-class]{GRanges} with slot
#' start holding the start position of the peak, slot end holding the end
#' position of the peak, slot space holding the chromosome location where the
#' peak is located, slot rownames holding the id of the peak. In addition, the
#' following variables are included.  \item{list("feature")}{id of the feature
#' such as ensembl gene ID} \item{list("insideFeature")}{upstream: peak resides
#' upstream of the feature; downstream: peak resides downstream of the feature;
#' inside: peak resides inside the feature; overlapStart: peak overlaps with
#' the start of the feature; overlapEnd: peak overlaps with the end of the
#' feature; includeFeature: peak include the feature entirely}
#' \item{list("distancetoFeature")}{distance to the nearest feature such as
#' transcription start site.  By default, the distance is calculated as the
#' distance between the start of the binding site and the TSS that is the gene
#' start for genes located on the forward strand and the gene end for genes
#' located on the reverse strand. The user can specify the location of peak and
#' location of feature for calculating this}
#' \item{list("start_position")}{start position of the feature such as gene}
#' \item{list("end_position")}{end position of the feature such as the gene}
#' \item{list("strand")}{1 or + for positive strand and -1 or - for negative
#' strand where the feature is located} \item{list("shortestDistance")}{The
#' shortest distance from either end of peak to either end the feature.  }
#' \item{list("fromOverlappingOrNearest")}{nearest: indicates this feature's
#' start (feature's end for features at minus strand) is closest to the peak
#' start; Overlapping: indicates this feature overlaps with this peak although
#' it is not the nearest feature start }
#' @author Lihua Julie Zhu, Jianhong Ou
#' @seealso \link{getAnnotation}, \link{findOverlappingPeaks},
#' \link{makeVennDiagram}, \link{addGeneIDs}, \link{peaksNearBDP},
#' \link{summarizePatternInPeaks}, \link{annoGR}, \link{annoPeaks}
#' @references 1. Zhu L.J. et al. (2010) ChIPpeakAnno: a Bioconductor package
#' to annotate ChIP-seq and ChIP-chip data. BMC Bioinformatics 2010,
#' 11:237doi:10.1186/1471-2105-11-237
#' 
#' 2. Zhu L (2013). "Integrative analysis of ChIP-chip and ChIP-seq dataset."
#' In Lee T and Luk ACS (eds.), Tilling Arrays, volume 1067, chapter 4, pp.
#' -19.  Humana Press. http://dx.doi.org/10.1007/978-1-62703-607-8_8
#' @keywords misc
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom BiocGenerics start end width strand
#' @importFrom S4Vectors mcols
#' @examples
#' 
#' 
#' #if (interactive()){
#'     ## example 1: annotate myPeakList by TxDb or EnsDb.
#'     data(myPeakList)
#'     library(ensembldb)
#'     library(EnsDb.Hsapiens.v75)
#'     annoData <- annoGR(EnsDb.Hsapiens.v75)
#'     annotatePeak = annotatePeakInBatch(myPeakList[1:6], AnnotationData=annoData)
#'     annotatePeak
#'     
#'     ## example 2: annotate myPeakList (GRanges) 
#'     ## with TSS.human.NCBI36 (Granges)
#'     data(TSS.human.NCBI36)
#'     annotatedPeak = annotatePeakInBatch(myPeakList[1:6], 
#'                                         AnnotationData=TSS.human.NCBI36)
#'     annotatedPeak
#'     
#'     ## example 3: you have a list of transcription factor biding sites from 
#'     ## literature and are interested in determining the extent of the overlap 
#'     ## to the list of peaks from your experiment. Prior calling the function 
#'     ## annotatePeakInBatch, need to represent both dataset as GRanges 
#'     ## where start is the start of the binding site, end is the end of the 
#'     ## binding site, names is the name of the binding site, space and strand 
#'     ## are the chromosome name and strand where the binding site is located.
#'     
#'     myexp <- GRanges(seqnames=c(6,6,6,6,5,4,4), 
#'                      IRanges(start=c(1543200,1557200,1563000,1569800,
#'                                      167889600,100,1000),
#'                              end=c(1555199,1560599,1565199,1573799,
#'                                    167893599,200,1200),
#'                              names=c("p1","p2","p3","p4","p5","p6", "p7")), 
#'                      strand="+")
#'     literature <- GRanges(seqnames=c(6,6,6,6,5,4,4), 
#'                           IRanges(start=c(1549800,1554400,1565000,1569400,
#'                                           167888600,120,800),
#'                                   end=c(1550599,1560799,1565399,1571199,
#'                                         167888999,140,1400),
#'                                   names=c("f1","f2","f3","f4","f5","f6","f7")),
#'                           strand=rep(c("+", "-"), c(5, 2)))
#'     annotatedPeak1 <- annotatePeakInBatch(myexp, 
#'                                           AnnotationData=literature)
#'     pie(table(annotatedPeak1$insideFeature))
#'     annotatedPeak1
#'     ### use toGRanges or rtracklayer::import to convert BED or GFF format
#'     ###  to GRanges before calling annotatePeakInBatch
#'     test.bed <- data.frame(space=c("4", "6"), 
#'                            start=c("100", "1000"),
#'                            end=c("200", "1100"), 
#'                            name=c("peak1", "peak2"))
#'     test.GR = toGRanges(test.bed)
#'     annotatePeakInBatch(test.GR, AnnotationData = literature)
#' #}
#' 
annotatePeakInBatch <-
    function (myPeakList, mart, featureType = c("TSS", "miRNA", "Exon"),
              AnnotationData, 
              output = c("nearestLocation", "overlapping", "both", 
                         "shortestDistance", "inside",
                         "upstream&inside", "inside&downstream",
                         "upstream", "downstream", 
                         "upstreamORdownstream", 
                         "nearestBiDirectionalPromoters"),
              multiple = c(TRUE,FALSE), maxgap = -1L,
              PeakLocForDistance=c("start","middle","end"),
              FeatureLocForDistance=c("TSS","middle","start","end", "geneEnd"),
              select=c("all", "first", "last", "arbitrary"),
              ignore.strand=TRUE, bindingRegion=NULL, ...)
    {
        if(output[1]=="nearestStart") output <- "nearestLocation"
        featureType = match.arg(featureType)
        PeakLocForDistance = match.arg(PeakLocForDistance)
        FeatureLocForDistance = match.arg(FeatureLocForDistance)
        output = match.arg(output)
        select = match.arg(select)
        multiple = multiple[1]
        
        if ((output == "overlapping" || output == "both")
            && select =="all" && multiple==FALSE) {
            warning("Please use select instead of multiple!")
            select = "first"
        }
        if(output=="upstream&inside"){
            if(FeatureLocForDistance!="TSS") {
                FeatureLocForDistance <- "TSS"
                warning("FeatureLocForDistance is set to TSS")
            }
            select <- "all"
        }
        if(output=="inside&downstream"){
            if(FeatureLocForDistance!="geneEnd"){
                FeatureLocForDistance <- "geneEnd"
                warning("FeatureLocForDistance is set to geneEnd")
            }
            select <- "all"
        }
        if (missing(myPeakList)) stop("Missing required argument myPeakList!")
        if (!is(myPeakList, "GRanges")) {
            stop("No valid myPeakList passed in. It needs to be GRanges object")
        }
        if (missing(AnnotationData)) {
            message("No AnnotationData as GRanges is passed in, 
                    so now querying biomart database for AnnotationData ....")
            if (missing(mart)) {
                stop("Error in querying biomart database. 
                     No valid mart object is passed in! 
                     Suggest call getAnnotation before calling 
                     annotatePeakInBatch")
            }
            if(!is(mart, "Mart")){
                if(is(mart, "GRanges")){
                    warning("AnnotationData is missing, and a GRanges is passed
                            to mart parameter. Trying to annotate peaks by 
                            the GRanges you input.")
                    AnnotationData <- mart
                    mart <- NULL
                }else{
                    stop("Error in querying biomart database. 
                     No valid mart object is passed in! 
                     Suggest call getAnnotation before calling 
                     annotatePeakInBatch")
                }
            }else{
                AnnotationData <- getAnnotation(mart, featureType = featureType)
                message("Done querying biomart database, start annotating ....
                        Better way would be calling getAnnotation before 
                        calling annotatePeakInBatch")
            }
        }
        if (!inherits(AnnotationData, 
                      c("GRanges", "annoGR"))) {
            stop("AnnotationData needs to be GRanges or annoGR object")
        }
        if(inherits(AnnotationData, "annoGR")){
            TSS.ordered <- as(AnnotationData, "GRanges")
        }else{
            TSS.ordered <- AnnotationData
        }
        nAnno <- length(TSS.ordered)
        
        rm(AnnotationData)
        rm(nAnno)
        
        if(length(bindingRegion)>1){
            message("Annotate peaks by annoPeaks, see ?annoPeaks for details.")
            ## FeatureLocForDistance=c("TSS","middle","start","end", "geneEnd")
            ##"startSite", "endSite", "fullRange"
            dots <- list(...)
            if(!"bindingType" %in% names(dots)){
                if(output %in% 
                   c("overlapping", 
                     "nearestBiDirectionalPromoters")){
                    bindingType <- switch(
                        output,
                        nearestBiDirectionalPromoters="nearestBiDirectionalPromoters",
                        overlapping=switch(FeatureLocForDistance,
                                           TSS="startSite",
                                           geneEnd="endSite",
                                           NA),
                        NA
                    )
                }
            }else{
              bindingType <- dots$bindingType
            }
            if(exists("bindingType") && !is.na(bindingType)){
                message("maxgap will be ignored.")
                dots$peak <- myPeakList
                dots$annoData <- TSS.ordered
                dots$bindingType <- bindingType
                dots$bindingRegion <- bindingRegion
                return(do.call(annoPeaks, dots))
            }else{
                message("bindingRegion will be ignored.")
            }
        }
        
        if (is.null(names(TSS.ordered))){
            names(TSS.ordered) <- formatC(1:length(TSS.ordered),
                                          width=nchar(length(TSS.ordered)),
                                          flag="0")
        }
        if (is.null(names(myPeakList))) {
            names(myPeakList) <- formatC(1:length(myPeakList),
                                         width = nchar(length(myPeakList)),
                                         flag = "0")
        }
        if(any(duplicated(names(myPeakList)))){
            warning("Found duplicated names in myPeakList. 
                    Changing the peak names ...")
            names(myPeakList) <- formatC(1:length(myPeakList),
                                         width = nchar(length(myPeakList)),
                                         flag = "0")
        }
        savedNames <- names(myPeakList)
        
        ##clear seqnames, the format should be chr+NUM
        ##TODO, how about the seqname not start with chr?
        TSS.ordered <- formatSeqnames(TSS.ordered)
        myPeakList <- formatSeqnames(myPeakList)
        if(!all(seqlevels(myPeakList) %in% seqlevels(TSS.ordered))){
            warning("not all the seqnames of myPeakList is 
                    in the AnnotationData.")
        }
        
        if(length(myPeakList)>10000){
            ##huge dataset
            myPeakList <- split(myPeakList, 
                                cut(1:length(myPeakList),
                                    ceiling(length(myPeakList)/5000)))
            myPeakList <- lapply(myPeakList, annotatePeakInBatch, 
                                 AnnotationData=TSS.ordered, 
                                 output = output, maxgap = maxgap,
                                 PeakLocForDistance=PeakLocForDistance,
                                 FeatureLocForDistance=FeatureLocForDistance,
                                 select=select,
                                 ignore.strand=ignore.strand)
            names(myPeakList) <- NULL
            myPeakList <- unlist(GRangesList(myPeakList))
            names(myPeakList) <- make.names(paste(myPeakList$peak, 
                                                  myPeakList$feature))
            ##myPeakList
            return(myPeakList)
        }
        ###STEP1 get nearst annotation for each peak, 
        ###use distanceToNearest(query, subject, ignore.strand=T/F, select)
        ## the distance got here is the shortestDistance
        ## select could only be arbitrary or all, 
        ## if it is "first" or "last", use "all" instead.
        ## if output=="nearest", annotation should only consider the the start point
        ##         ignore.strand <- all(strand(myPeakList)=="*") || 
        ##             all(strand(TSS.ordered)=="*") || 
        ##             all(strand(myPeakList)=="+")
        nsel <- ifelse(select %in% c("all", "first", "last"), 
                       "all", "arbitrary")
        featureGR <- TSS.ordered
        end(featureGR) <- 
            start(featureGR) <- 
            switch(FeatureLocForDistance,
                   TSS=ifelse(strand(featureGR)=="-", 
                              end(featureGR), 
                              start(featureGR)), 
                   geneEnd=ifelse(strand(featureGR)=="-", 
                                  start(featureGR), 
                                  end(featureGR)),
                   middle=round(rowMeans(cbind(start(featureGR), 
                                               end(featureGR)))),
                   start=start(featureGR),
                   end=end(featureGR) 
            )
        myPeaksGR <- myPeakList
        
        if(output=="nearestLocation"){
            dist <- as.data.frame(nearest(myPeaksGR, featureGR, 
                                          ignore.strand=ignore.strand, 
                                          select=nsel))            
            if(nrow(dist)==0) dist[1,] <- NA ##in case no match at all
            if(nsel=="arbitrary") {
                dist <- cbind(queryHits=seq_along(myPeakList), 
                              subjectHits=dist)
                colnames(dist) <- c("queryHits", "subjectHits")
            }
            dist$output <- rep("NearestLocation", nrow(dist))
        }
        if(output=="both"){
            distN <- as.data.frame(nearest(myPeaksGR, featureGR,
                                           ignore.strand=ignore.strand,
                                           select=nsel))
            if(nrow(distN)==0) distN[1,]<-NA
            if(nsel=="arbitrary") {
                distN <- cbind(queryHits=seq_along(myPeakList), 
                               subjectHits=distN)
                colnames(distN) <- c("queryHits", "subjectHits")
            }
            distO <- as.data.frame(findOverlaps(myPeakList, TSS.ordered,
                                                maxgap=maxgap,
                                                ignore.strand=ignore.strand,
                                                select=select))
            if(nrow(distO)==0) distO[1,]<-NA
            if(ncol(distO)==1){
                distO <- cbind(queryHits=seq_along(myPeakList), 
                               subjectHits=distO)
                colnames(distO) <- c("queryHits", "subjectHits")
            }
            distN$output <- rep("NearestLocation", nrow(distN))
            distN <- distN[!is.na(distN$subjectHits),]
            distO$output <- rep("Overlapping", nrow(distO))
            distO <- distO[!is.na(distO$subjectHits),]
            dist <- rbind(distN, distO)
            dist <- dist[!duplicated(dist[,c("queryHits", "subjectHits")]),
                         ,drop=FALSE]
            dist <- dist[order(dist$queryHits, dist$subjectHits),,drop=FALSE]
        }
        if(output=="overlapping"){
            dist <- as.data.frame(findOverlaps(myPeakList, TSS.ordered,
                                               maxgap=maxgap,
                                               ignore.strand=ignore.strand,
                                               select=select))
            if(nrow(dist)==0) dist[1,] <- NA
            if(ncol(dist)==1){
                dist <- cbind(queryHits=seq_along(myPeakList), subjectHits=dist)
                colnames(dist) <- c("queryHits", "subjectHits")
            }
            dist$output <- rep("Overlapping", nrow(dist))
        }
        if(output=="shortestDistance"){
            dist <- as.data.frame(nearest(myPeakList, TSS.ordered,
                                          ignore.strand=ignore.strand,
                                          select=nsel))
            if(nrow(dist)==0) dist[1,] <- NA ##in case no match at all
            if(nsel=="arbitrary") {
                dist <- cbind(queryHits=seq_along(myPeakList), subjectHits=dist)
                colnames(dist) <- c("queryHits", "subjectHits")
            }
            dist$output <- rep("shortestDistance", nrow(dist))
        }
        if(output=="upstream&inside"){
            #upstream
            featureGR <- TSS.ordered
            start <- ifelse(strand(featureGR)=="-", 
                            start(featureGR), 
                            start(featureGR)-max(maxgap, 1))
            width <- width(featureGR) + max(maxgap, 1)
            start(featureGR) <- start
            width(featureGR) <- width
            dist <- as.data.frame(findOverlaps(myPeakList, featureGR,
                                                ignore.strand=ignore.strand,
                                                select=select,
                                                type="any"))
            dist$output <- rep("Upstream&Inside", nrow(dist))
        }
        if(output=="upstream"){
            #upstream
            featureGR <- TSS.ordered
            start <- ifelse(strand(featureGR)=="-", 
                            end(featureGR)+1, 
                            start(featureGR)-max(maxgap, 1))
            width <- max(maxgap, 1)
            start(featureGR) <- start
            width(featureGR) <- width
            dist <- as.data.frame(findOverlaps(myPeakList, featureGR,
                                               ignore.strand=ignore.strand,
                                               select=select,
                                               type="any"))
            dist$output <- rep("Upstream", nrow(dist))
        }
        if(output=="inside&downstream"){
            #downstream
            featureGR <- TSS.ordered
            start <- ifelse(strand(featureGR)=="-", 
                            start(featureGR)-max(maxgap, 1), 
                            start(featureGR))
            width <- width(featureGR) + max(maxgap, 1)
            start(featureGR) <- start
            width(featureGR) <- width
            dist <- as.data.frame(findOverlaps(myPeakList, featureGR,
                                                ignore.strand=ignore.strand,
                                                select=select,
                                                type="any"))
            dist$output <- rep("Inside&Downstream", nrow(dist))
        }
        if(output=="downstream"){
            #downstream
            featureGR <- TSS.ordered
            start <- ifelse(strand(featureGR)=="-", 
                            start(featureGR)-max(maxgap, 1), 
                            end(featureGR)+1)
            width <- max(maxgap, 1)
            start(featureGR) <- start
            width(featureGR) <- width
            dist <- as.data.frame(findOverlaps(myPeakList, featureGR,
                                               ignore.strand=ignore.strand,
                                               select=select,
                                               type="any"))
            dist$output <- rep("Downstream", nrow(dist))
        }
        if(output=="upstreamORdownstream"){
            featureGR <- TSS.ordered
            start <- ifelse(strand(featureGR)=="-", 
                            end(featureGR)+1, 
                            start(featureGR)-max(maxgap, 1))
            width <- max(maxgap, 1)
            start(featureGR) <- start
            width(featureGR) <- width
            dist1 <- as.data.frame(findOverlaps(myPeakList, featureGR,
                                               ignore.strand=ignore.strand,
                                               select=select,
                                               type="any"))
            dist1$output <- rep("Upstream", nrow(dist1))
            featureGR <- TSS.ordered
            start <- ifelse(strand(featureGR)=="-", 
                            start(featureGR)-max(maxgap, 1), 
                            end(featureGR)+1)
            width <- max(maxgap, 1)
            start(featureGR) <- start
            width(featureGR) <- width
            dist2 <- as.data.frame(findOverlaps(myPeakList, featureGR,
                                               ignore.strand=ignore.strand,
                                               select=select,
                                               type="any"))
            dist2$output <- rep("Downstream", nrow(dist2))
            dist <- rbind(dist1, dist2)
        }
        if(output=="inside"){
            dist <- as.data.frame(findOverlaps(myPeakList, TSS.ordered,
                                               ignore.strand=ignore.strand,
                                               select=select,
                                               type="within"))
            dist$output <- rep("inside", nrow(dist))
        }
        if(output=="upstream2downstream"){
            featureGR <- TSS.ordered
            start(featureGR) <- 
                apply(cbind(start(featureGR) - maxgap, 1), 1, max)
            end(featureGR) <- end(featureGR) + maxgap
            dist <- as.data.frame(findOverlaps(myPeakList, featureGR,
                                               ignore.strand=ignore.strand,
                                               select=select,
                                               type="any"))
            dist$output <- rep("upstream2downstream", nrow(dist))
        }
        ##  nearest is NOT filtered by maxgap, is this should be changed?
        ##        dist <- dist[!is.na(dist$subjectHits), ]
        ##        distance <- distance(myPeakList[dist$queryHits],
        ##                             TSS.ordered[dist$subjectHits],
        ##                             ignore.strand=ignore.strand)
        ##        dist <- dist[abs(distance) <= maxgap, ]
        if(output=="nearestBiDirectionalPromoters" && length(bindingRegion)<=1){
          stop("If output is nearestBiDirectionalPromoters,",
               "please set bindingRegion. See ?annoPeaks for details.")
        }
        myPeakList.Hit <- 
            myPeakList[dist$queryHits[!is.na(dist$subjectHits)]]
        myPeakList.NA <- 
            myPeakList[!names(myPeakList) %in% names(myPeakList.Hit)]
        subjectHits <- 
            TSS.ordered[dist$subjectHits[!is.na(dist$subjectHits)]]
        mcols(subjectHits)$output <- 
            dist[!is.na(dist$subjectHits),"output"]
        #    myPeakList.Hit$distanceToNearest <- 
        #    dist$distance[!is.na(dist$subjectHits)]
        
        ###STEP2 get distance for each peak and nearest annotation by 
        ## distance(x, y). the distance is calculated by
        ##        PeakLocForDistance=c("start","middle","end"),
        ##        FeatureLocForDistance=c("TSS","middle","start",
        ##                              "end", "geneEnd")
        FeatureLoc <-
            switch(FeatureLocForDistance,
                   middle=as.integer(round(rowMeans(cbind(start(subjectHits), 
                                                          end(subjectHits))))),
                   start=start(subjectHits),
                   end=end(subjectHits),
                   geneEnd=as.integer(ifelse(strand(subjectHits)=="-", 
                                             start(subjectHits), 
                                             end(subjectHits))),
                   TSS=as.integer(ifelse(strand(subjectHits)=="-", 
                                         end(subjectHits), 
                                         start(subjectHits))))
        
        PeakLoc <- 
            switch(PeakLocForDistance,
                   start=start(myPeakList.Hit),
                   end=end(myPeakList.Hit),
                   middle=as.integer(round(rowMeans(
                       cbind(start(myPeakList.Hit), end(myPeakList.Hit))))))
        
        distancetoFeature <- as.numeric(ifelse(strand(subjectHits)=="-", 
                                               FeatureLoc - PeakLoc, 
                                               PeakLoc - FeatureLoc))
        
        ###STEP3 relationship between query and subject:
        ###   "inside", "overlapEnd", "overlapStart", 
        ###    "includeFeature", "upstream", "downstream"
        insideFeature <- getRelationship(myPeakList.Hit, subjectHits)
        myPeakList.Hit$peak <- names(myPeakList.Hit)
        myPeakList.Hit$feature <- names(subjectHits)
        myPeakList.Hit$start_position <- start(subjectHits)
        myPeakList.Hit$end_position <- end(subjectHits)
        myPeakList.Hit$feature_strand <- as.character(strand(subjectHits))
        
        myPeakList.Hit$insideFeature <- insideFeature[, "insideFeature"]
        myPeakList.Hit$distancetoFeature <- distancetoFeature
        
        myPeakList.Hit$shortestDistance <- 
            as.integer(insideFeature[,"shortestDistance"])
        
        myPeakList.Hit$fromOverlappingOrNearest <- subjectHits$output
        ##save oid for select == "first" or "last" filter
        myPeakList.Hit$oid <- insideFeature[, "ss"]
        ###combind data of NA with Hits
        if(length(myPeakList.NA)>0){
            myPeakList.NA$peak <- names(myPeakList.NA)
            for(ncol in c("feature",
                          "start_position",
                          "end_position",
                          "feature_strand",
                          "insideFeature",
                          "distancetoFeature",
                          "shortestDistance",
                          "fromOverlappingOrNearest",
                          "oid")){
                mcols(myPeakList.NA)[, ncol] <- NA
            }
            ## in case duplicated colnames
            colnames(mcols(myPeakList.NA)) <- colnames(mcols(myPeakList.Hit))
            myPeakList <- c(myPeakList.Hit, myPeakList.NA)
        }else{
            myPeakList <- myPeakList.Hit
        }
        ###output the results
        ##select=c("all", "first", "last", "arbitrary"))
        ##output = c("nearest", "overlapping", "both")
        ##if select %in% first or last, must filter the duplicate annotation
        ##if output="both", must filter the duplicate annotation
        ##if input is GRanges, output is GRanges
        ###the order of output should be same as input
        if(select=="first"){
            myPeakList <- 
                myPeakList[order(names(myPeakList), abs(myPeakList$oid))]
            myPeakList <- 
                myPeakList[!duplicated(names(myPeakList))]
        }
        if(select=="last"){
            myPeakList <- 
                myPeakList[order(names(myPeakList), -abs(myPeakList$oid))]
            myPeakList <- 
                myPeakList[!duplicated(names(myPeakList))]
        }
        ##remove column oid
        myPeakList$oid <- NULL
        ##re-order myPeakList as the original order
        oid <- 1:length(savedNames)
        names(oid) <- savedNames
        oid <- oid[names(myPeakList)]
        if(!any(is.na(oid))){
            myPeakList <- myPeakList[order(oid)]
        }
        

        if(output=="nearestLocation"){
            ## remove duplicate annotation, just keep the nearest one
            removeDuplicates <- function(gr){
                dup <- duplicated(gr$peak)
                if(any(dup)){
                    gr$oid <- 1:length(gr)
                    dup <- gr[dup]
                    gr.dup <- gr[gr$peak %in% dup$peak] 
                    ## bugs peak name must be different.
                    gr.NOTdup <- gr[!gr$peak %in% dup$peak]
                    gr.dup <- split(gr.dup, gr.dup$peak)
                    gr.dup <- lapply(gr.dup, function(.ele){
                        .ele[.ele$shortestDistance == 
                                 min(.ele$shortestDistance)]
                    })
                    gr.dup <- unlist(GRangesList(gr.dup))
                    gr <- c(gr.dup, gr.NOTdup)
                    gr <- gr[order(gr$oid)]
                    gr$oid <- NULL  
                }
                gr
            }
            myPeakList <- removeDuplicates(myPeakList)
        }

        
        names(myPeakList) <- 
            make.names(paste(myPeakList$peak, myPeakList$feature))

        ##myPeakList
        return(myPeakList)
    }
