binOverFeature <- function(..., annotationData=GRanges(),
                           select=c("all", "nearest"),
                           radius=5000L, nbins=50L,
                           minGeneLen=1L, aroundGene=FALSE, mbins=nbins, 
                           featureSite=c("FeatureStart", "FeatureEnd", 
                                         "bothEnd"),
                           PeakLocForDistance=c("all", "end", 
                                                "start", "middle"), 
                           FUN=sum, errFun=sd, xlab, ylab, main)
{
    ###check inputs
    PeaksList <- list(...)
    isGRangesList <- FALSE
    if(length(PeaksList)==1){
        if(is(PeaksList[[1]], "GRangesList")){
            PeaksList <- PeaksList[[1]]
            names <- names(PeaksList)
            isGRangesList <- TRUE
        }
    }
    ##save dots arguments names
    if(!isGRangesList){
        dots <- substitute(list(...))[-1]
        names <- unlist(sapply(dots, deparse))
    }
    
    n <- length(PeaksList)
    if(!n){
        stop("Missing required argument Peaks!")
    }else{
        nr <- ceiling(sqrt(n))
        nc <- ceiling(n/nr)
        op <- par(mfrow=c(nr, nc))
        on.exit(par(op))
    }
    if (missing(annotationData)) {
        stop("No AnnotationData as GRanges or annoGR is passed in.")
    }
    if(!inherits(annotationData, c("GRanges","RangedData", "annoGR")) || 
           length(annotationData)<1){
        stop("No AnnotationData as GRanges or annoGR is passed in.")
    }
    if(class(annotationData)=="annoGR")
        annotationData <- as(annotationData, "GRanges")
    annotationData <- unique(annotationData)
    if(class(annotationData)=="RangedData") 
        annotationData <- toGRanges(annotationData)
    if (!all(as.character(strand(annotationData)) %in% c("+", "-", "*")))
        stop("strands of annotationData must be +, - or *")
    select <- match.arg(select)
    featureSite <- match.arg(featureSite)
    PeakLocForDistance <- match.arg(PeakLocForDistance)
    if(!is.list(FUN)) FUN <- list(FUN)
    if(!is.list(errFun)) errFun <- list(errFun)
    lapply(FUN, function(fun){
        if(mode(fun)!="function")
            stop("The mode of FUN must be function. 
                 The FUN could be any function such as 
                 median, mean, sum, length, ...")
    })
    lapply(errFun, function(fun){
        if(mode(fun)!="function" & !is.numeric(fun))
            stop("The mode of errFun must be function. 
                 The errFun could be any function such as sd")
    })
    
    annotatedPeaksList<-lapply(PeaksList, function(Peaks){
        if (!inherits(Peaks, "GRanges")) {
            stop("No valid Peaks passed in. It needs to be GRanges object")
        }
        if(is.null(score(Peaks))){
            message("score of GRanges object is required for calculation. 
                    It will be the input of FUN. Setting score as 1.")
            Peaks$score <- 1
        }
        ## annotate the peaks
        if(select=="all"){
            annotatedPeaks <- 
                annotatePeakInBatch(Peaks, AnnotationData=annotationData,
                                    output="overlapping", maxgap = radius, 
                                    select="all") 
        }else{
            annotatedPeaks <- 
                annotatePeakInBatch(Peaks, AnnotationData=annotationData,
                                    output="nearestLocation", select="all")
        }
        ## filter the annotation
        annotatedPeaks <- annotatedPeaks[!is.na(annotatedPeaks$feature_strand)]
        annotatedPeaks <- 
            annotatedPeaks[
                as.numeric(as.character(annotatedPeaks$end_position))-
                    as.numeric(as.character(annotatedPeaks$start_position))+1>=
                    minGeneLen]
        
        ###if insideFeature==inside or includeFeature, 
        ###featureSite=="bothEnd.intergenic", the annotation should be removed
        #    if(featureSite=="bothEnd.intergenic") 
        #      annotatedPeaks <- 
        #        annotatedPeaks[!annotatedPeaks$insideFeature %in% 
        #                c("inside", "includeFeature"),]
        annotatedPeaks
    })
    
    plotErrBar <- function(x, y, err){
        yplus <- y+err
        yminus <- y-err
        segments(x, yminus, x, yplus)
        xcoord <- par()$usr[1:2]
        smidge <- 0.015*(xcoord[2]-xcoord[1])/2
        segments(x-smidge, yminus, x+smidge, yminus)
        segments(x-smidge, yplus, x+smidge, yplus)
    }
    
    if(missing(xlab)) {
        xlab <- if(aroundGene){
            paste("Bins from", featureSite)
        }else{
            paste("distance from", featureSite)
        }
    }
    if(missing(ylab)) ylab <- "Score" 
    if(missing(main)) main <- paste(names, "binding over", featureSite)
    
    binValue <- mapply(function(annotatedPeaks, fun, errfun,
                                xlab.ele, ylab.ele, main.ele){
        ##genelength
        genelen <- 
            annotatedPeaks$end_position - annotatedPeaks$start_position + 1
        ## change the function if it is length
        if(identical(fun, length)){
            fun <- sum
            annotatedPeaks$score <- 1
        }
        ###step 1 calculate the distance,
        strand <- annotatedPeaks$feature_strand=="-"
        if(PeakLocForDistance=="all"){
            ##dist, the start distance to feature loc, dist2, 
            ##the end distance to feature loc
            PeakLoc.start <- start(annotatedPeaks)
            PeakLoc.end <- end(annotatedPeaks)
            if(featureSite=="bothEnd"){
                ##bothEnd, TODO
                stop("Sorry, can not handle the combination of 
                     PeakLocForDistance=='all' AND featureSite=='bothEnd'")
            }else{
                FeatureLoc<-
                    switch(featureSite,
                           FeatureStart=ifelse(strand, 
                                               annotatedPeaks$end_position, 
                                               annotatedPeaks$start_position),
                           FeatureEnd=ifelse(strand, 
                                             annotatedPeaks$start_position, 
                                             annotatedPeaks$end_position),
                           0)
                dist1 <- ifelse(strand, 
                                FeatureLoc-PeakLoc.end, 
                                PeakLoc.start-FeatureLoc)
                dist2 <- ifelse(strand, 
                                FeatureLoc-PeakLoc.start, 
                                PeakLoc.end-FeatureLoc)
                score <- score(annotatedPeaks)
                weight <- (radius * mbins)/(genelen * nbins)
                if(aroundGene){
                    if(featureSite=="FeatureStart"){
                        ibin1 <- 
                            ifelse(dist1<0, 
                                   nbins+floor(dist1*nbins/radius),
                                   ifelse(dist1<genelen, 
                                          nbins+floor(dist1*mbins/genelen), 
                                          nbins+mbins+
                                              floor((dist1-genelen)*
                                                        nbins/radius)))
                        ibin2 <- 
                            ifelse(dist2<0, 
                                   nbins+floor(dist2*nbins/radius),
                                   ifelse(dist2<genelen, 
                                          nbins+floor(dist2*mbins/genelen), 
                                          nbins+mbins+
                                              floor((dist2-genelen)*
                                                        nbins/radius)))
                        b <- c(-1*c(nbins:1),0:(mbins+nbins-1))
                        distance2 <- -1*radius
                        distance1 <- genelen+radius
                    }else{##featureSite=="FeatureEnd"
                        ibin1 <- 
                            ifelse(dist1>=0, 
                                   nbins+mbins+floor(dist1*nbins/radius),
                                   ifelse(dist1>-1*genelen, 
                                          nbins+mbins+
                                              floor(dist1*mbins/genelen), 
                                          nbins+
                                              floor((dist1+genelen)*
                                                        nbins/radius)))
                        ibin2 <- 
                            ifelse(dist2>=0, 
                                   nbins+mbins+floor(dist2*nbins/radius),
                                   ifelse(dist2>-1*genelen, 
                                          nbins+mbins+
                                              floor(dist2*mbins/genelen), 
                                          nbins+
                                              floor((dist2+genelen)*
                                                        nbins/radius)))
                        b <- c(-1*c((mbins+nbins):1), 0:(nbins-1))
                        distance2 <- -1*(radius+genelen)
                        distance1 <- radius
                    }
                    numbins <- 2*nbins + mbins
                    b.type <- c(rep(FALSE, nbins), 
                                rep(TRUE, mbins), 
                                rep(FALSE, nbins))
                }else{
                    ibin1 <- round(nbins+floor(dist1*nbins/radius))
                    ibin2 <- round(nbins+floor(dist2*nbins/radius))
                    b <- c(-1*c(nbins:1), 0:(nbins-1))
                    numbins <- 2*nbins
                    distance2 <- -1*radius
                    distance1 <- radius
                    genelen <- 0
                    minGeneLen <- -1
                    b.type <- rep(FALSE, 2*nbins)
                }
                ids <- dist2>=distance2 & 
                    dist1<distance1 & 
                    genelen > minGeneLen
                score <- score[ids]
                weight <- weight[ids]
                ibin1 <- ibin1[ids]
                ibin2 <- ibin2[ids]
                ibin1[ibin1<0] <- 0
                ibin2[ibin2>numbins] <- numbins
                gps <- lapply(0:(numbins-1), function(.id){
                    ##split the scores
                    .idx <- ibin1<=.id & ibin2>=.id
                    if(b.type[.id+1]){
                        score[.idx] * weight[.idx]
                    }else{
                        score[.idx]
                    }
                })
                names(gps) <- formatC(1:numbins, 
                                      width=nchar(numbins), 
                                      flag="0")
            }
        }else{##PeakLocForDistance %in% start, middle, end
            PeakLoc <- 
                switch(PeakLocForDistance,
                       middle=round(rowMeans(cbind(start(annotatedPeaks), 
                                                   end(annotatedPeaks)))),
                       start=start(annotatedPeaks),
                       end=end(annotatedPeaks),
                       0)
            if(featureSite=="bothEnd"){##only consider outside of bothEnd
                FeatureStart=ifelse(strand, 
                                    annotatedPeaks$end_position, 
                                    annotatedPeaks$start_position)
                FeatureEnd=ifelse(strand, 
                                  annotatedPeaks$start_position, 
                                  annotatedPeaks$end_position)
                dist1 <- ifelse(strand, 
                                FeatureStart-PeakLoc, 
                                PeakLoc-FeatureStart) ##frome feature start
                dist2 <- ifelse(strand, 
                                FeatureEnd-PeakLoc, 
                                PeakLoc-FeatureEnd) ##from feature end
                score <- score(annotatedPeaks)
                if(aroundGene){
                    stop("Sorry, can not handle the combination of 
                         aroundGene==TRUE AND featureSite=='bothEnd'")
                }else{
                    ibin1 <- round(nbins+floor(dist1*nbins/radius))
                    ibin2 <- round(nbins+floor(dist2*nbins/radius))
                    score1 <- score[ibin1<nbins & ibin1>=0]
                    ibin1 <- ibin1[ibin1<nbins & ibin1>=0]
                    score2 <- score[ibin2>=nbins & ibin2<2*nbins]
                    ibin2 <- ibin2[ibin2>=nbins & ibin2<2*nbins]
                    numbins <- 2*nbins
                    b <- c(-1*c(nbins:1), 0:(nbins-1))
                    b.type <- rep(FALSE, length(b))
                    ##insert the empty bins
                    ibin1 <- formatC(ibin1, width=nchar(numbins), flag="0")
                    ibin2 <- formatC(ibin2, width=nchar(numbins), flag="0")
                    gps1 <- split(score1, ibin1)
                    gps2 <- split(score2, ibin2)
                    gps1 <- gps1[formatC(0:numbins, 
                                         width=nchar(numbins), 
                                         flag="0")]
                    gps2 <- gps2[formatC(0:numbins, 
                                         width=nchar(numbins), 
                                         flag="0")]
                    names(gps1) <- formatC(0:numbins, 
                                           width=nchar(numbins), 
                                           flag="0")
                    names(gps2) <- formatC(0:numbins, 
                                           width=nchar(numbins), 
                                           flag="0")
                    gps <- c(gps1[1:nbins], gps2[(nbins+1):(2*nbins)])
                }
            }else{##featureSite %in% FeatureStart, FeatureEnd
                FeatureLoc <-
                    switch(featureSite,
                           FeatureStart=ifelse(strand, 
                                               annotatedPeaks$end_position, 
                                               annotatedPeaks$start_position),
                           FeatureEnd=ifelse(strand, 
                                             annotatedPeaks$start_position, 
                                             annotatedPeaks$end_position),
                           0)
                dist1 <- ifelse(strand, 
                                FeatureLoc - PeakLoc, 
                                PeakLoc - FeatureLoc)
                weight <- (radius * mbins)/(genelen * nbins)
                if(aroundGene){
                    if(featureSite=="FeatureStart"){
                        ibin1 <- 
                            ifelse(dist1<0, 
                                   nbins+floor(dist1*nbins/radius),
                                   ifelse(dist1<genelen, 
                                          nbins+floor(dist1*mbins/genelen), 
                                          nbins+mbins+
                                              floor((dist1-genelen)*
                                                        nbins/radius)))
                        b <- c(-1*c(nbins:1),0:(mbins+nbins-1))
                        distance1 <- radius+genelen
                        distance2 <- -1*radius
                    }else{##featureSite=="FeatureEnd"
                        ibin1 <- 
                            ifelse(dist1>=0, 
                                   nbins+mbins+floor(dist1*nbins/radius),
                                   ifelse(dist1>-1*genelen, 
                                          nbins+mbins+
                                              floor(dist1*mbins/genelen), 
                                          nbins+
                                              floor((dist1+genelen)*
                                                        nbins/radius)))
                        b <- c(-1*c((mbins+nbins):1), 0:(nbins-1))
                        distance1 <- radius
                        distance2 <- -1*(radius+genelen)
                    }
                    numbins <- 2*nbins + mbins
                    b.type <- c(rep(FALSE, nbins), 
                                rep(TRUE, mbins), 
                                rep(FALSE, nbins))
                }else{
                    ibin1 <- round(nbins+floor(dist1*nbins/radius))
                    numbins <- 2*nbins
                    b <- c(-1*c(nbins:1), 0:(nbins-1))
                    genelen <- 0
                    minGeneLen <- -1
                    distance1 <- radius
                    distance2 <- -1*radius
                    b.type <- rep(FALSE, 2 * nbins)
                }
                ids <- dist1>=distance2 & 
                    dist1<distance1 & 
                    genelen > minGeneLen
                score <- score(annotatedPeaks)[ids]
                weight <- weight[ids]
                ibin1 <- ibin1[ids]
                ibin1[ibin1<0] <- 0
                ibin1[ibin1>numbins] <- numbins
                score[b.type[ibin1+1]] <- 
                    score[b.type[ibin1+1]] * weight[b.type[ibin1+1]]
                ibin1 <- formatC(ibin1, 
                                 width=nchar(numbins), 
                                 flag="0")
                gps <- split(score, ibin1)
                gps <- gps[formatC(0:numbins, 
                                   width=nchar(numbins),
                                   flag="0")]
                names(gps) <- formatC(0:numbins, 
                                      width=nchar(numbins), 
                                      flag="0")
                ##insert the empty bins
                gps <- gps[1:numbins]
            }
        }
        gps <- lapply(gps, function(.ele){
            if(is.null(.ele[1])){
                    0
                }else{
                    .ele
                }
            })
        value <- unlist(lapply(gps, fun))
        std <- if(mode(errfun)=="function") unlist(lapply(gps, errfun)) else rep(errfun, length(gps))
        std[is.na(std)] <- 0
        ##plot the figure
        ylim.min <- min(value[!is.na(value)] - std[!is.na(value)])
        ylim.max <- max(value[!is.na(value)] + std[!is.na(value)])
        ylim.dis <- (ylim.max - ylim.min)/20
        blabel <- if(aroundGene){
            c(seq.int(-radius, -radius/nbins, length.out=nbins)+radius/nbins/2,
              1:mbins,
              seq.int(0, radius-radius/nbins, length.out=nbins)+radius/nbins/2)
        }else{
            seq.int(-radius, radius-radius/nbins, length.out=2*nbins) +
                radius/nbins/2
        }
        if(aroundGene){
            plot(b, value, 
                 ylim=c(ylim.min-ylim.dis, ylim.max+ylim.dis),
                 xlab=xlab.ele, 
                 ylab=ylab.ele, 
                 main=main.ele,
                 xaxt="n")
            b.type.at <- b[which(b.type)]
            b.type.at <- 
                b.type.at[c(1, length(b.type.at))] + c(-.5, .5)
            abline(v=b.type.at, lty=2)
            b.at <- c(b[1]-.5, b[nbins]+.5, b[nbins+mbins]+.5, b[length(b)]+.5)
            b.label <- c(-1 * radius, "Feature Start", "Feature End", radius)
            axis(1, at=b.at, labels=b.label)
            if(!all(std==0)) plotErrBar(b, value, std)
        }else{
            plot(blabel, value, 
                 ylim=c(ylim.min-ylim.dis, ylim.max+ylim.dis),
                 xlab=xlab.ele, 
                 ylab=ylab.ele, 
                 main=main.ele)
            if(!all(std==0)) plotErrBar(blabel, value, std)
        }
        
        names(value) <- blabel
        value
    }, annotatedPeaksList, FUN, errFun, xlab, ylab, main)
    
    colnames(binValue) <- names
    ###output statistics
    return(invisible(binValue))
}
