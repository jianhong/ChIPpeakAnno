overLapWithHg38Tss <- function(input.bw.path,output.file.dir,dd) {
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  data(TSS.human.GRCh38)
  
  FeatureLocForDistance=c("TSS","middle","start","end", "geneEnd")
  
  FeatureLocForDistance = "TSS"
  
  TSS.ordered <- TSS.human.GRCh38
  
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
  
  dd = 5000
  
  r.start <- start(featureGR)-dd
  r.end <-  start(featureGR)+dd
  
  featureGR.1 <- featureGR
  
  start(featureGR.1) <- r.start
  end(featureGR.1) <- r.end
  
  zz.temp <- as.data.frame(featureGR.1)
  
  #zz.temp <- as.data.frame(TSS.human.GRCh38)
  
  chrOrder<-c(1:22,"X","Y")
  
  zz.temp$seqnames <- factor(zz.temp$seqnames, levels=chrOrder)
  zz.temp.1 <- zz.temp[order(zz.temp$seqnames,zz.temp$start),]
  zz.temp.2 <- zz.temp.1[-which(is.na(zz.temp.1$seqnames)),]
  
  z <- data.frame(paste("chr",zz.temp.2[,1],sep=""),zz.temp.2[,c(2:3,5)],row.names(zz.temp.2))
  colnames(z)=c("seqnames","start","end","strand","label")
  z <- GRanges(z)
  
  z.L <- list(Hg38Tss=z)
  
  class_pattern <- "Hg38Tss"
  data <- z.L
  
  tss <- names(data)[grep(class_pattern,names(data),ignore.case = T)]
  print(tss)
  
  tss.bed <- data[which(names(data) %in% tss)]
  
  # dd <- 5000
  #   
  #   u <- 1
  #   features <- tss.bed[[u]]
  #   
  #   wid <- width(features)
  #   
  #   feature.center <- features
  #   feature.recentered <- feature.center
  #   
  #   start(feature.center) <- start(features) + floor(wid/2)
  #   
  #   width(feature.center) <- 1
  #   
  #   start(feature.recentered) <- start(feature.center) - dd
  #   end(feature.recentered) <- end(feature.center) + dd
  
  feature.recentered <- tss.bed[[1]] 
  
  input.bw.path <- "~/Aimin/DropboxUmass/Aimin/Project/Magnolia/PositiveControl"
  
  files <- dir(input.bw.path, "bigWig")
  if(.Platform$OS.type != "windows"){
    cvglists <- sapply(file.path(input.bw.path,files), import, 
                       format="BigWig", 
                       which=feature.recentered, 
                       as="RleList")
  }else{## rtracklayer can not import bigWig files on Windows
    load(file.path(input.bw.path, "cvglist.rds"))
  }
  
  names(cvglists) <- gsub(".bw", "", files)
  
  feature.center <- reCenterPeaks(tss.bed[[1]], width=1)
  
  sig <- featureAlignedSignal(cvglists, feature.center, 
                              upstream=5000, downstream=5000)
  
  names(sig) <- "positive_control_GSE95979"
  
  #keep <- rowSums(sig[[1]]) > 0
  #sig <- sig[[1]][keep, ]
  #feature.center <- feature.center[keep]
  
  heatmap <- ChIPpeakAnno::featureAlignedHeatmap(sig, feature.center, 
                                                 upstream=5000, downstream=5000)
  
  
  orderMat4HeatMap <- function(Mat4Ht,scaling=F,orderBy) {
    
    Mat4Ht <- Mat4Ht[-which(apply(Mat4Ht,1,sd)==0),]
    
    if(scaling==T)
      {y <- t(scale(t(Mat4Ht), center=TRUE, scale=TRUE))}
    else
    {y <- Mat4Ht}
    
    density.sum <- apply(y,1,orderBy)
    sig.sorted<- y[order(density.sum,decreasing = T),]        
    z.mat.0 <- sig.sorted
    z.mat.0
  }
  
  z.mat.0 <- orderMat4HeatMap(sig[[1]],scaling=T,orderBy=max)
  
  z.mat.0 <- orderMat4HeatMap(sig[[1]],scaling=T,orderBy=sum)
  
  print(rownames(z.mat.0))
  
  cols <- colorpanel(100,"white","lightblue","blue")
  
  ht <- Heatmap(z.mat.0, name = "coverage",
                #col = rev(redblue(30))[-seq(35, 35)],  
                #col = rev(redblue(30)),
                # col = colorpanel(30,"white","white","red"),
                col = cols,
                show_row_name = FALSE,
                cluster_columns = F,
                cluster_rows = F,column_names_gp = gpar(fontsize = 6))
  
  ht
  
  
  
  featureAlignedDistribution(sig, feature.center, 
                             upstream=5000, downstream=5000,
                             type="l")
  
  sig.normalized <- lapply(sig, function(u){
    y <- t(scale(t(u), center=TRUE, scale=TRUE))
    y
  })
  
  heatmap <- ChIPpeakAnno::featureAlignedHeatmap(sig.normalized, feature.center, 
                                                 upstream=5000, downstream=5000)
  
  #keep <- rowSums(sig[[2]]) > 0
  #sig <- sapply(sig, function(.ele) .ele[keep, ], simplify = FALSE)
  #feature.center <- feature.center[keep]
  
  # heatmap <- featureAlignedHeatmap(sig, feature.center, 
  #                                  upstream=5000, downstream=5000,
  #                                  lower.extreme=c(-400.8430,-612.4750,-456.1310,-690.1795,
  #                                                  -475.8967,-582.0910),
  #                                  upper.extreme=c(860.7876,777.8088,596.1618,
  #                                                  291.7805,253.2713,238.3586))
  # 
  library(gplots)
  cols <- colorpanel(100,"white","lightblue","blue")
  heatmap <- featureAlignedHeatmap(sig.normalized, feature.center, 
                                   upstream=5000, color=cols,downstream=5000)
  
  featureAlignedDistribution(sig.normalized, feature.center, 
                             upstream=5000, downstream=5000,
                             type="l")
  
  input.bw.path <- "/Users/aiminyan/Aimin/DropboxUmass/Aimin/Project/Magnolia/HSePeaks"
  
  files <- dir(input.bw.path, "bw")
  
  files.1 <- files[c(1,3,4,6,8,10,12,13,14,16)]
  
  files.2 <- files.1[c(1,3,8,9)]  
  
  TssPlot <- function(.Platform, input.bw.path, files.2, feature.center, feature.recentered) {
    if(.Platform$OS.type != "windows"){
      cvglists <- sapply(file.path(input.bw.path,files.2), import, 
                         format="BigWig", 
                         which=feature.recentered, 
                         as="RleList")
    }else{## rtracklayer can not import bigWig files on Windows
      load(file.path(input.bw.path, "cvglist.rds"))
    }
    
    names(cvglists) <- gsub(".bw", "", files.2)
    
    sig <- featureAlignedSignal(cvglists, feature.center, 
                                upstream=5000, downstream=5000)
    
    names(sig) <- gsub("1218HS_","",names(sig))
    
    featureAlignedDistribution(sig, feature.center, 
                               upstream=5000, downstream=5000,
                               type="l")
    
    heatmap <- ChIPpeakAnno::featureAlignedHeatmap(sig, feature.center, 
                                                   upstream=5000, downstream=5000)
  }
  
  TssPlot(.Platform, input.bw.path, files.2, feature.center, feature.recentered)
  
  files.1 <- files[c(1,3,4,6,8,10,12,13,14,16)]
  
  files.3 <- files.1[c(2,4,5,6,7,10)]  
  
  TssPlot(.Platform, input.bw.path, files.3, feature.center, feature.recentered)
  
  files.4 <- files.1[c(6)]
  TssPlot(.Platform, input.bw.path, files.4, feature.center, feature.recentered)
  
  files.4 <- files.1[c(6)]
  TssPlot(.Platform, input.bw.path, files.4, feature.center, feature.recentered)
  
  files.4 <- files.1[c(1)]
  TssPlot(.Platform, input.bw.path, files.4, feature.center, feature.recentered)
  
  sig.normalized <- lapply(sig, function(u){
    y <- t(scale(t(u), center=TRUE, scale=TRUE))
    y
  })
  
  cols <- colorpanel(100,"white","lightblue","blue")
  heatmap <- featureAlignedHeatmap(sig.normalized, feature.center, 
                                   upstream=5000, color=cols,downstream=5000)
  
  featureAlignedDistribution(sig.normalized, feature.center, 
                             upstream=5000, downstream=5000,
                             type="l")
  
  featureAlignedDistribution(sig, feature.center, 
                             upstream=5000, downstream=5000,
                             type="l")
  
  cvglists <- list(A=RleList(chr1=Rle(sample.int(5000, 100), 
                                      sample.int(300, 100))), 
                   B=RleList(chr1=Rle(sample.int(5000, 100), 
                                      sample.int(300, 100))))
  feature.gr <- GRanges("chr1", IRanges(seq(1, 4900, 100), width=100))
  featureAlignedDistribution(cvglists, feature.gr, zeroAt=50, type="l")
  
}
