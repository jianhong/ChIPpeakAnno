annotatePeakInBatch <-
    function (myPeakList, mart, featureType = c("TSS", "miRNA", "Exon"),
              AnnotationData, 
              output = c("nearestLocation", "overlapping", "both", 
                         "shortestDistance", "inside",
                         "upstream&inside", "inside&downstream",
                         "upstream", "downstream", 
                         "upstreamORdownstream"),
              multiple = c(TRUE,FALSE), maxgap = 0L,
              PeakLocForDistance=c("start","middle","end"),
              FeatureLocForDistance=c("TSS","middle","start","end", "geneEnd"),
              select=c("all", "first", "last", "arbitrary"),
              ignore.strand=TRUE)
    {
        if(output[1]=="nearestStart") output <- "nearestLocation"
        featureType = match.arg(featureType)
        PeakLocForDistance = match.arg(PeakLocForDistance)
        FeatureLocForDistance = match.arg(FeatureLocForDistance)
        output = match.arg(output)
        select = match.arg(select)
        
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
        if (!inherits(myPeakList, c("RangedData", "GRanges"))) {
            stop("No valid myPeakList passed in. It needs to be RangedData or GRanges object")
        }
        if(inherits(myPeakList, "RangedData")){
            myPeakList <- RangedData2GRanges(myPeakList)
            flagRD <- TRUE
        }else{
            flagRD <- FALSE
        }
        if (missing(AnnotationData)) {
            message("No AnnotationData as RangedData or GRanges is passed in, so now querying biomart database for AnnotationData ....")
            if (missing(mart) || class(mart) != "Mart") {
                stop("Error in querying biomart database. No valid mart object is passed in! Suggest call getAnnotation before calling annotatePeakInBatch")
            }
            AnnotationData <- getAnnotation(mart, featureType = featureType)
            message("Done querying biomart database, start annotating ....Better way would be calling getAnnotation before calling annotatePeakInBatch")
        }
        if (!inherits(AnnotationData, c("RangedData", "GRanges"))) {
            stop("AnnotationData needs to be RangedData or GRanges object")
        }
        if (inherits(AnnotationData, "RangedData")) {
            TSS.ordered <- RangedData2GRanges(AnnotationData)
        }else{
            TSS.ordered <- AnnotationData
        }
        nAnno <- length(TSS.ordered)
        #TSS.ordered <- unique(TSS.ordered)
        #if(length(TSS.ordered)!=nAnno) message("AnnotationData has duplicated ranges.")
        rm(AnnotationData)
        rm(nAnno)
        
        ### rownames of RangedData changed into names of GRanges
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
            warning("Found duplicated names in myPeakList. Changing the peak names ...")
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
            warning("not all the seqnames of myPeakList is in the AnnotationData.")
        }
        
        if(length(myPeakList)>10000){
            ##huge dataset
            myPeakList <- split(myPeakList, cut(1:length(myPeakList), ceiling(length(myPeakList)/5000)))
            myPeakList <- lapply(myPeakList, annotatePeakInBatch, 
                                 AnnotationData=TSS.ordered, 
                                 output = output, maxgap = maxgap,
                                 PeakLocForDistance=PeakLocForDistance,
                                 FeatureLocForDistance=FeatureLocForDistance,
                                 select=select,
                                 ignore.strand=ignore.strand)
            names(myPeakList) <- NULL
            myPeakList <- unlist(GRangesList(myPeakList))
            if(flagRD){
                names(myPeakList) <- paste(myPeakList$peak, myPeakList$feature)
                ##output RangedData
                myPeakList <- as(myPeakList, "RangedData")
                ##myPeakList <- myPeakList[!is.na(myPeakList$feature),]
            }else{
                names(myPeakList) <- make.names(paste(myPeakList$peak, myPeakList$feature))
            }
            ##myPeakList
            return(myPeakList)
        }
        ###STEP1 get nearst annotation for each peak, use distanceToNearest(query, subject, ignore.strand=T/F, select)
        ## the distance got here is the shortestDistance
        ## select could only be arbitrary or all, if it is "first" or "last", use "all" instead.
        ## if output=="nearest", annotation should only consider the the start point
#         ignore.strand <- all(strand(myPeakList)=="*") || 
#             all(strand(TSS.ordered)=="*") || 
#             all(strand(myPeakList)=="+")
        nsel <- ifelse(select %in% c("all", "first", "last"), "all", "arbitrary")
        featureGR <- TSS.ordered
        end(featureGR) <- start(featureGR) <- switch(FeatureLocForDistance,
                                                     TSS=ifelse(strand(featureGR)=="-", end(featureGR), start(featureGR)), 
                                                     geneEnd=ifelse(strand(featureGR)=="-", start(featureGR), end(featureGR)),
                                                     middle=round(rowMeans(cbind(start(featureGR), end(featureGR)))),
                                                     start=start(featureGR),
                                                     end=end(featureGR) 
            )
        myPeaksGR <- myPeakList
#         end(myPeaksGR) <- start(myPeaksGR) <-switch(PeakLocForDistance,
#                                                     start=ifelse(strand(myPeakList)=="-",
#                                                                  end(myPeakList),
#                                                                  start(myPeakList)),
#                                                     end=ifelse(strand(myPeakList)=="-",
#                                                                start(myPeakList),
#                                                                end(myPeakList)),
#                                                     middle=round(rowMeans(cbind(start(myPeakList), end(myPeakList)))))
        if(output=="nearestLocation"){
            dist <- as.data.frame(nearest(myPeaksGR, featureGR, ignore.strand=ignore.strand, select=nsel))            
            if(nrow(dist)==0) dist[1,] <- NA ##in case no match at all
            if(nsel=="arbitrary") {
                dist <- cbind(queryHits=1:length(myPeakList), subjectHits=dist)
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
                distN <- cbind(queryHits=1:length(myPeakList), subjectHits=distN)
                colnames(distN) <- c("queryHits", "subjectHits")
            }
            distO <- as.data.frame(findOverlaps(myPeakList, TSS.ordered,
                                                maxgap=maxgap,
                                                ignore.strand=ignore.strand,
                                                select=select))
            if(nrow(distO)==0) distO[1,]<-NA
            if(ncol(distO)==1){
                distO <- cbind(queryHits=1:length(myPeakList), subjectHits=distO)
                colnames(distO) <- c("queryHits", "subjectHits")
            }
            distN$output <- rep("NearestLocation", nrow(distN))
            distN <- distN[!is.na(distN$subjectHits),]
            distO$output <- rep("Overlapping", nrow(distO))
            distO <- distO[!is.na(distO$subjectHits),]
            dist <- rbind(distN, distO)
            dist <- dist[!duplicated(dist[,c("queryHits", "subjectHits")]),,drop=FALSE]
            dist <- dist[order(dist$queryHits, dist$subjectHits),,drop=FALSE]
        }
        if(output=="overlapping"){
            dist <- as.data.frame(findOverlaps(myPeakList, TSS.ordered,
                                               maxgap=maxgap,
                                               ignore.strand=ignore.strand,
                                               select=select))
            if(nrow(dist)==0) dist[1,] <- NA
            if(ncol(dist)==1){
                dist <- cbind(queryHits=1:length(myPeakList), subjectHits=dist)
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
                dist <- cbind(queryHits=1:length(myPeakList), subjectHits=dist)
                colnames(dist) <- c("queryHits", "subjectHits")
            }
            dist$output <- rep("shortestDistance", nrow(dist))
        }
        if(output=="upstream&inside"){
            #upstream
            featureGR <- TSS.ordered
            start <- ifelse(strand(featureGR)=="-", start(featureGR), start(featureGR)-max(maxgap, 1))
            width <- width(featureGR) + max(maxgap, 1)
            start(featureGR) <- start
            width(featureGR) <- width
            dist <- as.data.frame(findOverlaps(myPeakList, featureGR,
                                                maxgap=0,
                                                ignore.strand=ignore.strand,
                                                select=select,
                                                type="any"))
            dist$output <- rep("Upstream&Inside", nrow(dist))
        }
        if(output=="upstream"){
            #upstream
            featureGR <- TSS.ordered
            start <- ifelse(strand(featureGR)=="-", end(featureGR)+1, start(featureGR)-max(maxgap, 1))
            width <- max(maxgap, 1)
            start(featureGR) <- start
            width(featureGR) <- width
            dist <- as.data.frame(findOverlaps(myPeakList, featureGR,
                                               maxgap=0,
                                               ignore.strand=ignore.strand,
                                               select=select,
                                               type="any"))
            dist$output <- rep("Upstream", nrow(dist))
        }
        if(output=="inside&downstream"){
            #downstream
            featureGR <- TSS.ordered
            start <- ifelse(strand(featureGR)=="-", start(featureGR)-max(maxgap, 1), start(featureGR))
            width <- width(featureGR) + max(maxgap, 1)
            start(featureGR) <- start
            width(featureGR) <- width
            dist <- as.data.frame(findOverlaps(myPeakList, featureGR,
                                                maxgap=0,
                                                ignore.strand=ignore.strand,
                                                select=select,
                                                type="any"))
            dist$output <- rep("Inside&Downstream", nrow(dist))
        }
        if(output=="downstream"){
            #downstream
            featureGR <- TSS.ordered
            start <- ifelse(strand(featureGR)=="-", start(featureGR)-max(maxgap, 1), end(featureGR)+1)
            width <- max(maxgap, 1)
            start(featureGR) <- start
            width(featureGR) <- width
            dist <- as.data.frame(findOverlaps(myPeakList, featureGR,
                                               maxgap=0,
                                               ignore.strand=ignore.strand,
                                               select=select,
                                               type="any"))
            dist$output <- rep("Downstream", nrow(dist))
        }
        if(output=="upstreamORdownstream"){
            featureGR <- TSS.ordered
            start <- ifelse(strand(featureGR)=="-", end(featureGR)+1, start(featureGR)-max(maxgap, 1))
            width <- max(maxgap, 1)
            start(featureGR) <- start
            width(featureGR) <- width
            dist1 <- as.data.frame(findOverlaps(myPeakList, featureGR,
                                               maxgap=0,
                                               ignore.strand=ignore.strand,
                                               select=select,
                                               type="any"))
            dist1$output <- rep("Upstream", nrow(dist1))
            featureGR <- TSS.ordered
            start <- ifelse(strand(featureGR)=="-", start(featureGR)-max(maxgap, 1), end(featureGR)+1)
            width <- max(maxgap, 1)
            start(featureGR) <- start
            width(featureGR) <- width
            dist2 <- as.data.frame(findOverlaps(myPeakList, featureGR,
                                               maxgap=0,
                                               ignore.strand=ignore.strand,
                                               select=select,
                                               type="any"))
            dist2$output <- rep("Downstream", nrow(dist2))
            dist <- rbind(dist1, dist2)
        }
        if(output=="inside"){
            dist <- as.data.frame(findOverlaps(myPeakList, TSS.ordered,
                                               maxgap=0,
                                               ignore.strand=ignore.strand,
                                               select=select,
                                               type="within"))
            dist$output <- rep("inside", nrow(dist))
        }
        if(output=="upstream2downstream"){
            featureGR <- TSS.ordered
            start(featureGR) <- apply(cbind(start(featureGR) - maxgap, 1), 1, max)
            end(featureGR) <- end(featureGR) + maxgap
            dist <- as.data.frame(findOverlaps(myPeakList, featureGR,
                                               maxgap=0,
                                               ignore.strand=ignore.strand,
                                               select=select,
                                               type="any"))
            dist$output <- rep("upstream2downstream", nrow(dist))
        }
##        nearest is NOT filtered by maxgap, is this should be changed?
##        dist <- dist[!is.na(dist$subjectHits), ]
##        distance <- distance(myPeakList[dist$queryHits],
##                             TSS.ordered[dist$subjectHits],
##                             ignore.strand=ignore.strand)
##        dist <- dist[abs(distance) <= maxgap, ]
        myPeakList.Hit <- myPeakList[dist$queryHits[!is.na(dist$subjectHits)]]
        myPeakList.NA <- myPeakList[!names(myPeakList) %in% names(myPeakList.Hit)]
        subjectHits <- TSS.ordered[dist$subjectHits[!is.na(dist$subjectHits)]]
        subjectHits$output <- dist[!is.na(dist$subjectHits),"output"]
        #    myPeakList.Hit$distanceToNearest <- dist$distance[!is.na(dist$subjectHits)]
        
        ###STEP2 get distance for each peak and nearest annotation by distance(x, y)
        ## the distance is calculated by
        ##        PeakLocForDistance=c("start","middle","end"),
        ##        FeatureLocForDistance=c("TSS","middle","start","end", "geneEnd")
        FeatureLoc<-switch(FeatureLocForDistance,
                           middle=as.integer(round(rowMeans(cbind(start(subjectHits), end(subjectHits))))),
                           start=start(subjectHits),
                           end=end(subjectHits),
                           geneEnd=as.integer(ifelse(strand(subjectHits)=="-", start(subjectHits), end(subjectHits))),
                           TSS=as.integer(ifelse(strand(subjectHits)=="-", end(subjectHits), start(subjectHits))))
        
        PeakLoc<-switch(PeakLocForDistance,
                        start=start(myPeakList.Hit),
                        end=end(myPeakList.Hit),
                        middle=as.integer(round(rowMeans(cbind(start(myPeakList.Hit), end(myPeakList.Hit))))))
        
        distancetoFeature <- as.numeric(ifelse(strand(subjectHits)=="-", FeatureLoc - PeakLoc, PeakLoc - FeatureLoc))
        
        ###STEP3 relationship between query and subject:
        ###   "inside", "overlapEnd", "overlapStart", "includeFeature", "upstream", "downstream"
        insideFeature <- getRelationship(myPeakList.Hit, subjectHits)
        myPeakList.Hit$peak <- names(myPeakList.Hit)
        myPeakList.Hit$feature <- names(subjectHits)
        myPeakList.Hit$start_position <- start(subjectHits)
        myPeakList.Hit$end_position <- end(subjectHits)
        myPeakList.Hit$feature_strand <- as.character(strand(subjectHits))
        
        myPeakList.Hit$insideFeature <- insideFeature[, "insideFeature"]
        myPeakList.Hit$distancetoFeature <- distancetoFeature
        
        myPeakList.Hit$shortestDistance <- as.integer(insideFeature[,"shortestDistance"])
        
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
            myPeakList <- c(myPeakList.Hit, myPeakList.NA)
        }else{
            myPeakList <- myPeakList.Hit
        }
        ###output the results
        ##select=c("all", "first", "last", "arbitrary"))
        ##output = c("nearest", "overlapping", "both")
        ##if select %in% first or last, must filter the duplicate annotation
        ##if output="both", must filter the duplicate annotation
        ###output should be a RangedData or GRanges?
        ##if input is RangedData, output is RangedData
        ##if input is GRanges, output is GRanges
        ###the order of output should be same as input
        if(select=="first"){
            myPeakList <- myPeakList[order(names(myPeakList), abs(myPeakList$oid))]
            myPeakList <- myPeakList[!duplicated(names(myPeakList))]
        }
        if(select=="last"){
            myPeakList <- myPeakList[order(names(myPeakList), -abs(myPeakList$oid))]
            myPeakList <- myPeakList[!duplicated(names(myPeakList))]
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
        

        if(output=="nearestLocation"){## remove duplicate annotation, just keep the nearest one
            removeDuplicates <- function(gr){
                dup <- duplicated(gr$peak)
                if(any(dup)){
                    gr$oid <- 1:length(gr)
                    dup <- gr[dup]
                    gr.dup <- gr[gr$peak %in% dup$peak] ## bugs peak name must be different.
                    gr.NOTdup <- gr[!gr$peak %in% dup$peak]
                    gr.dup <- split(gr.dup, gr.dup$peak)
                    gr.dup <- lapply(gr.dup, function(.ele){
                        .ele[.ele$shortestDistance == min(.ele$shortestDistance)]
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

        if(flagRD){
            names(myPeakList) <- paste(myPeakList$peak, myPeakList$feature)
            ##output RangedData
            myPeakList <- as(myPeakList, "RangedData")
            ##myPeakList <- myPeakList[!is.na(myPeakList$feature),]
        }else{
            names(myPeakList) <- make.names(paste(myPeakList$peak, myPeakList$feature))
        }
        ##myPeakList
        return(myPeakList)
    }
