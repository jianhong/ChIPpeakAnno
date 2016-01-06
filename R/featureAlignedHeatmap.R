featureAlignedHeatmap <- 
    function(cvglists, feature.gr, upstream, downstream, 
             zeroAt, n.tile=100,
             annoMcols=c(), sortBy=names(cvglists)[1],
             color=colorRampPalette(c("yellow", "red"))(50),
             lower.extreme, upper.extreme,
             margin=c(0.1, 0.01, 0.15, 0.1), gap=0.01, 
             newpage=TRUE, gp=gpar(fontsize=10),
             ...){
    if(!missing(lower.extreme)){
        stopifnot(is.numeric(lower.extreme))
    }
    if(!missing(upper.extreme)){
        stopifnot(is.numeric(upper.extreme))
    }
    stopifnot(class(gp)=="gpar")
    stopifnot(is.logical(newpage))
    stopifnot(class(feature.gr)=="GRanges")
    
    grWidr <- unique(width(feature.gr))
    if(missing(upstream) || missing(downstream)){
        if(length(grWidr)!=1){
            stop("The width of feature.gr is not identical.")
        }
        if(missing(zeroAt)) {
            zeroAt <- 0.5
            message("zero is set as the center of the feature.gr")
        }
    }else{
        if(!is.numeric(upstream) || !is.numeric(downstream)){
            stop("upstream and downstream must be integers")
        }
        if(upstream<0 || downstream<0){
            stop("upstream and downstream must be not less than 0")
        }
        upstream <- as.integer(upstream)
        downstream <- as.integer(downstream)
        if(length(grWidr)!=1){
            start(feature.gr) <- start(feature.gr)+floor(grWidr/2)
            width(feature.gr) <- 1
            warning("feature.gr is set to the center of feature.gr")
        }
        zeroAt <- upstream/(upstream + downstream)
        end(feature.gr) <- start(feature.gr) + downstream
        start(feature.gr) <- start(feature.gr) - upstream
        grWidr <- unique(width(feature.gr))
    }
    stopifnot(is.numeric(zeroAt))
    stopifnot(zeroAt>=0)
    if(zeroAt<=1){
        zero <- round(grWidr*zeroAt)
    }else{
        zero <- round(zeroAt)
    }
    
    grWid <- c(0, grWidr) - zero
    grWidLab <- grid.pretty(grWid)
    grWidAt <- (grWidLab+zero)/grWidr
    if(inherits(cvglists, c("SimpleRleList", "RleList", "CompressedRleList"))){
        cvglistsName <- substitute(deparse(cvglists))
        cvglists <- list(cvglists)
        names(cvglists) <- cvglistsName
    }
    if(class(cvglists)!="list"){
        stop("cvglists must be output of featureAlignedSignal or 
             a list of SimpleRleList or RleList")
    }
    cls <- sapply(cvglists, class)
    if(all(cls=="matrix")){
        cov <- cvglists
        if(ncol(cov[[1]])!=n.tile){
            stop("n.tile must keep same as featureAligendSignal.")
        }
    }else{
        if(any(!cls %in% c("SimpleRleList", "RleList", "CompressedRleList")))
            stop("cvglists must be a list of SimpleRleList or RleList")
        cov <- featureAlignedSignal(cvglists, feature.gr, n.tile=n.tile)
    }
    
    annoMcols <- annoMcols[annoMcols %in% colnames(mcols(feature.gr))]
    if(length(sortBy)>0){
        covTotalCoverage <- 
            do.call(cbind, lapply(cov, rowSums, na.rm=TRUE))
        colnames(covTotalCoverage) <- names(cvglists)
        if(length(annoMcols)>0){
            sortMatrix <- mcols(feature.gr)[, annoMcols, drop=FALSE]
            for(i in 1:ncol(sortMatrix)){
                sortMatrix[,i] <- 
                    factor(as.character(sortMatrix[,i]), 
                           levels=rev(unique(as.character(sortMatrix[, i]))))
            }
            sortMatrix <- 
                cbind(sortMatrix,
                      data.frame(covTotalCoverage[, sortBy, drop=FALSE]))
        }else{
            sortMatrix <- data.frame(covTotalCoverage[, sortBy, drop=FALSE])
        }
        names(sortMatrix) <- paste("V", 1:ncol(sortMatrix))
        soid <- do.call(order, c(as.list(sortMatrix), decreasing = TRUE))
        cov <- lapply(cov, function(.ele) .ele[soid, ])
        feature.gr <- feature.gr[soid]
    }
    
    #covert to color
    if(missing(lower.extreme)){
        lower.extreme <- rep(0, length(cvglists))
    }
    if(missing(upper.extreme)){
        lim=list()
    }else{
        upper.extreme <- 
            rep(upper.extreme, length(cvglists))[1:length(cvglists)]
        lim=as.list(data.frame(rbind(lower.extreme, upper.extreme)))
        names(lim) <- NULL
    }
    
    if(length(lim)==0){
        lim <- sapply(cov, range, na.rm=TRUE, simplify=FALSE)
    }
    if(length(lim)==1){
        lim <- rep(lim, length(cov))
        legend <- 1
    }else{
        legend <- length(cov)
    }
    if(any(sapply(lim, length)<2)){
        stop("each limit should be a numeric vector with length 2")
    }
    cov <- mapply(function(.ele, .lim){
        .lim <- seq(min(.lim), max(.lim), length.out=length(color)+1)
        .lim[1] <- min(.lim[1], 0)
        .lim[length(.lim)] <- max(max(.lim), .Machine$integer.max)
        t(apply(.ele, 1, 
                cut, breaks=.lim, labels=color, include.lowest=TRUE))
    }, cov, lim, SIMPLIFY=FALSE)
    # library(grid)
    # grid.raster
    if(newpage) grid.newpage()
    allGrob <- NULL
    height=1-margin[3]-margin[1]
    width=1-margin[4]-margin[2]
    ## draw y labels
    areColors <- function(x) {
        sapply(x, function(X) {
            tryCatch(is.matrix(col2rgb(X)), 
                     error = function(e) FALSE)
        })
    }
    colorGroup <- 
        list(c("#295AA5", "#52BDBD", "#FFAD6B", "#F784A5", "#CEB55A", 
               "#31B5D6", "#BD2119", "#F79463", "#FFE608", "#D69419", 
               "#29B5CE", "#737BB5"), 
             c("#4A4221", "#A5BD42", "#8CC584", "#CEE6D6", "#C5D64A", 
               "#FFF79C", "#84AD31", "#DBDECE", "#D6AD4A", "#4A4229", 
               "#637B7B", "#525229", "#F7F7AD", "#EFAD63"), 
             c("#6B5252", "#DE6B73", "#000000", "#4A4221", "#3AA55A", 
               "#DEE6C5", "#102119", "#DEE6CE", "#738CC5", "#31317B", 
               "#F7F7AD", "#8C9CA5", "#732163", "#001000"), 
             c("#192919", "#8C946B", "#6B1921", "#F7F7AD", "#4A4229", 
               "#CEA54A", "#000000", "#6B736B", "#524A29", "#423119", 
               "#213129", "#7B4229", "#7B2121", "#C57B4A"), 
             c("#9D2932", "#EFDEB0", "#789262", "#494166", "#A88462", 
               "#D9B612", "#177CB0", "#F3F8F1", "#424B50", "#B36D61", 
               "#C3272B", "#A98175", "#D4F2E8", "#FF3300", "#76664D", 
               "#ED5736", "#A3E2C5", "#415065", "#D7ECF1"), 
             c("#009C8C", "#0094AD", "#BDDEC5", "#00ADCE", "#DEFF19", 
               "#D6EFB5", "#EFF7E6", "#42CE73", "#42A5CE", "#BDFFF7", 
               "#DEFFDE", "#4A6BCE"), 
             c("#D6B5BD", "#EFD6BD", "#FFE6D6", "#DEE6CE", "#F7F7C5", 
               "#D6E6B5", "#73BD8C", "#BDDEDE", "#DEBD8C", "#E6A5BD", 
               "#F7ADAD", "#FFE6D6", "#F7F7D6", "#E6BDBD"), 
             c("#E6E6BD", "#637342", "#425A52", "#298473", "#BDDEC5", 
               "#527384", "#73847B", "#CED6AD", "#94D6CE", "#527384", 
               "#4A6B5A", "#318442", "#C5DEDE"), 
             c("#F7846B", "#F79C9C", "#CEE6E6", "#FFEFBD", "#EFEFBD", 
               "#FFE6DE", "#FFBDA5", "#FFEF00", "#63BD84", "#D6E6AD", 
               "#FFF7D6", "#FFFFEF", "#FFE6DE", "#FFBD9C", "#EF9C94"), 
             c("#4A4229", "#212108", "#6B8484", "#EFE6BD", "#E68C3A", 
               "#A56352", "#6B2921", "#525229", "#844229", "#8C524A", 
               "#E68C3A", "#D6AD4A", "#212108", "#524A3A"), 
             c("#312184", "#AD2919", "#EF3129", "#D69C4A", "#FF6B00", 
               "#6B4A31", "#CE3121", "#AD2919", "#5A3A5A", "#292931", 
               "#3A3184", "#9C5229", "#EFCE19"), 
             c("#6373B5", "#FF5200", "#4AF54A", "#FFCE10", "#EF3A21", 
               "#EF4A94", "#5A63AD", "#FFEDA5", "#EF3121", "#FFC510", 
               "#FF5200"), 
             c("#FFB6C1", "#DC143C", "#DB7093", "#FF69B4", "#FF1493", 
               "#C71585", "#DA70D6", "#FF00FF", "#8B008B", "#9400D3", 
               "#4B0082", "#FFA07A", "#FF7F50", "#FF4500", "#E9967A", 
               "#FF6347", "#F08080", "#CD5C5C", "#FF0000", "#A52A2A", 
               "#B22222", "#8B0000"), 
             c("#7B68EE", "#6A5ACD", "#483D8B", "#E6E6FA", "#0000FF", 
               "#191970", "#4169E1", "#6495ED", "#B0C4DE", "#1E90FF", 
               "#87CEFA", "#00BFFF", "#00FFFF"), 
             c("#008B8B", "#48D1CC", "#7FFFD4", "#66CDAA", "#00FA9A", 
               "#00FF7F", "#3CB371", "#2E8B57", "#90EE90", "#32CD32", 
               "#00FF00", "#008000", "#006400", "#ADFF2F", "#556B2F", 
               "#9ACD32", "#6B8E23"), 
             c("#FAFAD2", "#FFFF00", "#808000", "#BDB76B", "#FFFACD", 
               "#EEE8AA", "#F0E68C", "#FFD700", "#DAA520", "#B8860B", 
               "#F5DEB3", "#FFA500", "#FFEBCD", "#FFDEAD", "#D2B48C", 
               "#FF8C00", "#CD853F", "#8B4513"),
             c("#F5F5F5", "#DCDCDC", "#D3D3D3", "#C0C0C0", "#A9A9A9", 
               "#808080", "#696969", "#5C5C5C", "#4D4D4D", "#333333", 
               "#1A1A1A", "#000000"))
    if(length(annoMcols)>0){
        ci <- 1
        ytop <- 0.75-margin[3]
        mcols <- mcols(feature.gr)
        vp.annoMcols <- viewport(x=.5, y=height/2+margin[1], 
                                 width=1, height=height,
                                 name="vp.annoMcols")
        for(i in annoMcols){
            mc <- as.character(mcols[, i])
            mc.name <- ifelse(is.numeric(i), colnames(mcols)[i], i)
            if(!all(areColors(mc))){
                mc.label <- unique(mc)
                mc <- factor(mc, levels=mc.label)
                mc.label.color <- rep(colorGroup[[ci]], 
                                      length(mc.label))[1:length(mc.label)]
                levels(mc) <- mc.label.color
                mc <- as.character(mc)
                ci <- ci+1
                ## mc is color
                ## mc.label is legends
                ## plot legned
                ht.annoMcols.legend <- length(mc.label) * 0.03
                wd.annoMcols.legend <- margin[4]
                vp.annoMcols.legend <- 
                    viewport(x=1-wd.annoMcols.legend+gap, 
                             y=ytop-ht.annoMcols.legend/2, 
                             width=wd.annoMcols.legend, 
                             height=ht.annoMcols.legend,
                             just=0,
                             name="vp.annoMcols.legend")
                annoMcols.legend <- legendGrob(labels=mc.label, 
                                               ncol=1,
                                               pch=15,
                                               gp=gpar(col=mc.label.color))
                allGrob <- gList(allGrob, 
                                 gTree(children=gList(annoMcols.legend),
                                       vp=vp.annoMcols.legend,
                                       name=paste("gTree.annoMcols.legend",
                                                  mc.name,
                                                  sep=".")))
                ytop <- ytop-ht.annoMcols.legend-gap
            }
            
            vp.annoMcols.sub <- viewport(x=width - .0075,
                                         y=.5, 
                                         width=.015,
                                         height=1,
                                         name=paste("vp.annoMcols",
                                                    mc.name,
                                                    sep="."))
            raster.annoMcols.sub <- 
                rasterGrob(mc,
                           width=1, height=1,
                           interpolate=FALSE,
                           name=paste("raster.annoMcols",
                                      mc.name,
                                      sep="."),
                           vp=vp.annoMcols.sub)
            xaxis.annoMcols.sub <- 
                xaxisGrob(at=.5,
                          label=mc.name,
                          gp=gp,
                          edits=gEdit("labels", rot=90),
                          name=paste("xaxis.annoMcols",
                                     mc.name,
                                     sep="."),
                          vp=vp.annoMcols.sub)
            allGrob <- gList(allGrob, 
                             gTree(children=gList(raster.annoMcols.sub,
                                                  xaxis.annoMcols.sub),
                                   vp=vp.annoMcols,
                                   name=paste("gTree.annoMcols",
                                              mc.name,
                                              sep=".")))
            width <- width - .015
        }
        width <- width - gap
    }
    
    vp.heatmap <- viewport(x=width/2+margin[2], y=height/2+margin[1], 
                           width=width, height=height,
                           name="vp.heatmap")
    l <- length(cov)
    wid <- 1/l
    x <- -wid/2
    y <- 1
    for(i in 1:l){
        x <- x+wid
        x.name <- names(cvglists)[i]
        vp.heatmap.sub <- viewport(x=x, width=wid-gap, height=y,
                                   name=paste("vp.heatmap",
                                              x.name,
                                              sep="."))
        raster.heatmap.sub <- rasterGrob(cov[[i]], 
                                         width=1, height=1, 
                                         interpolate=FALSE,
                                         name=paste("raster.heatmap",
                                                    x.name,
                                                    sep="."),
                                         vp=vp.heatmap.sub)
        
        grWidLab.this <- grWidLab
        if(i==1){
            if(i!=l){
                if(length(grWidAt)>1){
                    grWidLab.this[length(grWidLab)] <- " "
                }
            }
        }else{
            if(i==l){##i!=1
                if(length(grWidAt)>1){
                    grWidLab.this[1] <- " "
                }
            }else{## middle
                if(length(grWidAt)>2){
                    grWidLab.this[1] <- " "
                    grWidLab.this[length(grWidLab)] <- " "
                }
            }
        }
        xaxis.heatmap.sub <- xaxisGrob(at=grWidAt,
                                       label=grWidLab.this,
                                       gp=gp,
                                       edits=gEdit("labels", rot=90),
                                       name=paste("xaxis.heatmap",
                                                  x.name,
                                                  sep="."),
                                       vp=vp.heatmap.sub)
        allGrob <- gList(allGrob, 
                         gTree(children=gList(raster.heatmap.sub,
                                              xaxis.heatmap.sub),
                               vp=vp.heatmap,
                               name=paste("gTree.heatmap",
                                          x.name,
                                          sep=".")))
    }
    
    ## add heatmap legend
    if(legend==1){
        vp.legend <- viewport(x=1.01-margin[4]+gap, y=.875-margin[3], 
                              width=.02, height=.25,
                              name="vp.legend")
        raster.legend <- rasterGrob(matrix(rev(color), ncol=1), 
                                    width=1, height=1,
                                    name="raster.legend")
        .lim <- lim[[1]]
        label <- grid.pretty(.lim)
        ll <- length(label)
        if(ll>=2){
            if(ll%%2==0){
                label <- label[c(1, ll)]
            }else{
                label <- label[c(1, ceiling(ll/2), ll)]
            }
            at <- (label-min(.lim))/abs(diff(.lim))
            yaxis.legend <- yaxisGrob(at=at,
                                      label=label,
                                      main=FALSE,
                                      gp=gp,
                                      name="yaxis.legend")
        }
        allGrob <- gList(allGrob, 
                         gTree(children=gList(raster.legend,
                                              yaxis.legend),
                               vp=vp.legend,
                               name="gTree.legend"))
        
        ## draw x labels
        vp.xlabels <- viewport(x=width/2+margin[2], 
                               y=1-margin[3]/2+gap, 
                               width=width, height=margin[3]-gap,
                               name="vp.xlabels")
        x <- -wid/2
        h <- .5
        for(i in 1:l){
            x <- x+wid
            x.name <- names(cvglists)[i]
            text.xlabels.sub <- 
                textGrob(label = x.name,
                         gp=gp,
                         just=1, hjust=0,
                         vp=viewport(x=x, y=.5, width=wid, height=1, 
                                     name=paste("vp.xlabels", x.name, sep=".")),
                         name=paste("text.xlabels", x.name, sep="."))
            allGrob <- gList(allGrob, 
                             gTree(children=gList(text.xlabels.sub),
                                   vp=vp.xlabels,
                                   name=paste("gTree.xlabels",
                                              x.name,
                                              sep=".")))
        }
    }else{
        vp.legend <- viewport(x=width/2+margin[2], 
                              y=1-margin[3]/3+gap, 
                              width=width, height=margin[3]*2/3-gap,
                              name="vp.legend")
        
        vp.xlabels <- viewport(x=width/2+margin[2], 
                               y=1-margin[3]*5/6+gap, 
                               width=width, height=margin[3]/3-gap,
                               name="vp.xlabels")
        
        x <- -wid/2
        h <- .5
        for(i in 1:l){
            x <- x+wid
            x.name=names(cvglists)[i]
            vp.legend.sub <- 
                viewport(x=x, y=0.25, width=wid/2, height=0.2,
                         name=paste("vp.legend",
                                    x.name,
                                    sep="."))
            raster.legend.sub <- 
                rasterGrob(matrix(color, nrow=1), 
                           width=1, height=1,
                           name=paste("raster.legend", x.name, sep="."),
                           vp=vp.legend.sub)
            .lim <- lim[[i]]
            label <- grid.pretty(.lim)
            ll <- length(label)
            if(ll>=2){
                if(ll%%2==0){
                    label <- label[c(1, ll)]
                }else{
                    label <- label[c(1, ceiling(ll/2), ll)]
                }
                at <- (label-min(.lim))/abs(diff(.lim))
                xaxis.legned.sub <- xaxisGrob(at=at,
                                              label=label,
                                              main=FALSE,
                                              gp=gp,
                                              name=paste("xaxis.legend",
                                                         x.name,
                                                         sep="."),
                                              vp=vp.legend.sub)
            }
            allGrob <- gList(allGrob, 
                             gTree(children=gList(raster.legend.sub,
                                                  xaxis.legned.sub),
                                   vp=vp.legend,
                                   name=paste("gTree.legend",
                                              x.name,
                                              sep=".")))
            
            ## draw x labels
            text.xlabels.sub <- 
                textGrob(label = x.name,
                         gp=gp,
                         vp=viewport(x=x, y=.5, width=wid, height=1, 
                                     name=paste("vp.xlabels", x.name, sep=".")),
                         name=paste("text.xlabels", x.name, sep="."))
            allGrob <- gList(allGrob, 
                             gTree(children=gList(text.xlabels.sub),
                                   vp=vp.xlabels,
                                   name=paste("gTree.xlabels",
                                              x.name,
                                              sep=".")))
            
        }
    }    
    
    ## consider how to return all and could be redraw by grid.draw.
    grid.draw(allGrob)
    return(invisible(allGrob))
}