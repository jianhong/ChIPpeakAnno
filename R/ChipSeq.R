#' AnnotatePeakUMASS
#'
#' @param input.file.dir 
#' @param input.file.pattern 
#' @param output.file.dir 
#' @param genome 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' input.file.dir='/Users/aiminyan/Aimin/DropboxUmass/NADfinder/BedFiles'
#' 
#' input.file.dir='~/Aimin/DropboxUmass/NADfinder/BedFiles/Aizhan\ F121-9\ comparisons'
#' 
#' input.file.pattern="*.bed$"
#' 
#' Output:
#' output.file.dir="/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output"
#' output.file.dir="/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output/OverLapPeakAdded"
#' output.file.dir="/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output/OverLapPeakAdded2"
#' output.file.dir="/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results"
#' 
#' output.file.dir="/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_9_5_2018"
#' 
#' output.file.dir="/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_9_6_2018/XL"
#' output.file.dir="/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_9_6_2018/nonXL"
#' output.file.dir="/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_9_10_2018"
#' 
#' output.file.dir="/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_9_10_2018/F121_9-vs-XL_MEF"
#' 
#' output.file.dir="/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_9_27_2018" 
#' 
#' AnnotatePeakUMASS(input.file.dir,input.file.pattern,output.file.dir,genome="Mm")
#' 
AnnotatePeakUMASS <- function(input.file.dir,input.file.pattern,output.file.dir,genome) {
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  file.name.4 <- list.files(input.file.dir,pattern=input.file.pattern,all.files = TRUE,full.names = TRUE,recursive = FALSE,include.dirs = TRUE)
  
  re.out<-lapply(file.name.4,function(u){
    
    if(length(grep("Peric",u)) ==1){
      peaks=read.table(u)
    }else if(length(grep("GSM",u)) ==1){
      peaks=read.table(u)
    }else
    {
      peaks=read.table(u,skip=1)
    }
    colnames(peaks)[1:3]= c("chr","start","end")
    peaks=toGRanges(peaks)
    peaks
  })
  
  names(re.out) <- gsub(" ","_",tools::file_path_sans_ext(basename(file.name.4)))
  
  pat1 <- "ciLAD"
  pat2 <- "nonXL_MEF"
  pat3 <- "XL_MEF"
  
  index.pat1 <- grep(paste("^",pat1, sep=""), names(re.out))
  index.pat2 <- grep(paste("^",pat2, sep=""), names(re.out))
  index.pat3 <- grep(paste("^",pat3, sep=""), names(re.out))
  
  ol <- findOverlapsOfPeaks(re.out[c(1,3)])
 
  makeNewListAndVenn <- function(re.out, index.pat1, index.pat2, output.file.dir, pat1, pat2) {
    
    pdf(file.path(output.file.dir,paste0(pat1,"-vs-",pat2,"_peak_overlap.pdf")))
    ol.by.min <- makeVennDiagram(re.out[c(index.pat1,index.pat2)], NameOfPeaks=c(pat1,pat2),
                                 totalTest=10000,scaled=F, euler.d=F,fill = c("red","blue"),
                                 alpha = 0.50,
                                 label.col = c(rep("black",3)),
                                 cex = 2,
                                 fontfamily = "serif",
                                 fontface = "bold",
                                 cat.col = c("red","blue"),connectedPeaks = "min")
    dev.off()
    
    ol.ciLAD.nonXL <- findOverlapsOfPeaks(re.out[c(index.pat1,index.pat2)])
    
    addOverLapPeak2List <- function(ol.ciLAD.nonXL, re.out) {
      overlapName <- gsub("///","-overlap-",names(ol.ciLAD.nonXL$peaklist)[3])
      re.out[[overlapName]] <- ol.ciLAD.nonXL$peaklist[[3]]
      re.out
    }
    
    re.out <- addOverLapPeak2List(ol.ciLAD.nonXL, re.out)
    re.out
  }
  
  ol.by.min <- makeVennDiagram(re.out[c(1,3)], NameOfPeaks=c("ciLAD","MEF_LAD"),
                               totalTest=10000,scaled=F, euler.d=F,fill = c("red","blue"),
                               alpha = 0.50,
                               label.col = c(rep("black",3)),
                               cex = 2,
                               fontfamily = "serif",
                               fontface = "bold",
                               cat.col = c("red","blue"),connectedPeaks = "min")
  
  ol.ciLAD.XL.MEF_LAD <- findOverlapsOfPeaks(re.out[c(1,4,3)])
  
  orderPeakAndOutPut <- function(peakAll,output.file.dir,output.file.name,outHeader= FALSE) {
    chrOrder<-c(paste("chr",1:19,sep=""),"chrX","chrY")
    
    peakAll$seqnames <- factor(peakAll$seqnames, levels=chrOrder)
    
    peakAll2 <- peakAll[order(peakAll$seqnames,peakAll$start),]
    
    if(outHeader == TRUE){
    write.table(peakAll2,file = file.path(output.file.dir,paste0(output.file.name,".txt")),append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T)
    }else{
    write.table(peakAll2,file = file.path(output.file.dir,paste0(output.file.name,".bed")),append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = F)
    }
  }
  
  overlap.with.ciLAD.but.XL_MEF_only_peaks <- ol.ciLAD.XL.MEF_LAD$overlappingPeaks[[3]][,c(11:13)]
  overlap.with.ciLAD.but.XL_MEF_only_peaks <- unique(overlap.with.ciLAD.but.XL_MEF_only_peaks)
  overlap.with.ciLAD.but.XL_MEF_only_peaks
  output.file.name ="overlap.with.ciLAD.but.XL_MEF_only_peaks"
  orderPeakAndOutPut(overlap.with.ciLAD.but.XL_MEF_only_peaks,output.file.dir,output.file.name)
    
  getOverLap <- function(ol.ciLAD.XL.MEF_LAD,index_overlap,output.file.dir) {
    X_name <- names(ol.ciLAD.XL.MEF_LAD$overlappingPeaks)
    
    overlappedPeaks <- data.frame(ol.ciLAD.XL.MEF_LAD$overlappingPeaks[[index_overlap]])
    names(overlappedPeaks)
    
    XX_name <- gsub("///",",", X_name[index_overlap])
  
    XXX_name <- unlist(strsplit(XX_name,","))
    
    colnames(overlappedPeaks)[which(colnames(overlappedPeaks)=="peaks1")] <- XXX_name[1]
    colnames(overlappedPeaks)[which(colnames(overlappedPeaks)=="peaks2")] <- XXX_name[2]
    
    overlap.peak.region <- apply(overlappedPeaks, 1, function(u){
      
      seqnames <- u["seqnames"]
      start <- max(u["start"],u["start.1"])
      end <- min(u["end"],u["end.1"])
      
      x <- data.frame(seqnames.olp=seqnames,start.olp=start,end.olp=end)      
      x
      
    })
    overlap.peak.region <- do.call(rbind.data.frame,overlap.peak.region)
    
    overlappedPeaks <- data.frame(overlappedPeaks,overlap.peak.region)
    
    overlappedPeaks <- data.frame(overlappedPeaks,overlappedBasePairs=as.numeric(as.character(overlappedPeaks$end.olp))
                                  -as.numeric(as.character(overlappedPeaks$start.olp))+1)
    
   
    output.file.name <- paste0(XXX_name[1],"-overlap-",XXX_name[2])
    
    output.file.name <- gsub("<","_", output.file.name)
    
    orderPeakAndOutPut(overlappedPeaks,output.file.dir,output.file.name,outHeader=T)
    overlappedPeaks
  }
  
  ciLAD.XL <- getOverLap(ol.ciLAD.XL.MEF_LAD,3,output.file.dir)
  getOverLap(ol.ciLAD.XL.MEF_LAD,2,output.file.dir)
  getOverLap(ol.ciLAD.XL.MEF_LAD,1,output.file.dir)
  
  ciLAD.part.in.ciLAD.XL <- unique(toGRanges(ciLAD.XL[,c(2,3,4)]))
  
  re.out[["ciLAD.part.in.ciLAD.XL"]] <- ciLAD.part.in.ciLAD.XL
  
  ol.ciLAD_XL_MEFLAD <- findOverlapsOfPeaks(re.out[c(6,3)])
  
  overlap.with.MEF_LAD.but.non_XL_MEF_only_peaks <- ol.ciLAD_XL_MEFLAD$overlappingPeaks[[1]][,c(2:4)]
  overlap.with.MEF_LAD.but.non_XL_MEF_only_peaks <- unique(overlap.with.MEF_LAD.but.non_XL_MEF_only_peaks)
  dim(overlap.with.MEF_LAD.but.non_XL_MEF_only_peaks)
  peak.in.ciLAD <- ol.ciLAD_XL_MEFLAD$uniquePeaks[grep("ciLAD.part.in.ciLAD.XL",names(ol.ciLAD_XL_MEFLAD$uniquePeaks)),]
  length(peak.in.ciLAD)
  
  temp <- ciLAD.XL[,c(11,12,13)]
  colnames(temp) <- c("seqnames","start","end")
  XL.part.in.ciLAD.XL <- unique(toGRanges(temp))
  
  re.out[["XL.part.in.ciLAD.XL"]] <- XL.part.in.ciLAD.XL
  
  ol.ciLAD_XL_MEFLAD <- findOverlapsOfPeaks(re.out[c(5,3)])
  
  overlap.with.MEF_LAD.but.non_XL_MEF_only_peaks <- ol.ciLAD_XL_MEFLAD$overlappingPeaks[[1]][,c(2:4)]
  overlap.with.MEF_LAD.but.non_XL_MEF_only_peaks <- unique(overlap.with.MEF_LAD.but.non_XL_MEF_only_peaks)
  dim(overlap.with.MEF_LAD.but.non_XL_MEF_only_peaks)
  
  peak.in.ciLAD <- ol.ciLAD_XL_MEFLAD$uniquePeaks[grep("ciLAD.XL",names(ol.ciLAD_XL_MEFLAD$uniquePeaks)),]
  length(peak.in.ciLAD)
  
  
  overlap.with.MEF_LAD.but.non_XL_MEF_only_peaks <- ol.ciLAD_XL_MEFLAD$overlappingPeaks[[1]][,c(8:10)]
  overlap.with.MEF_LAD.but.non_XL_MEF_only_peaks <- unique(overlap.with.MEF_LAD.but.non_XL_MEF_only_peaks)
  dim(overlap.with.MEF_LAD.but.non_XL_MEF_only_peaks)
  
  
  peak.in.ciLAD <- ol.ciLAD_XL_MEFLAD$uniquePeaks[grep("XL.part.in.ciLAD.XL",names(ol.ciLAD_XL_MEFLAD$uniquePeaks)),]
  length(peak.in.ciLAD)
  
  peak.in.ciLAD <- unique(ol.ciLAD_XL_MEFLAD$uniquePeaks[grep("Peric-Hupkes-2010-MEF-LADs-mm10",names(ol.ciLAD_XL_MEFLAD$uniquePeaks)),])
  length(peak.in.ciLAD)
  
  
  
  
  overlap.with.MEF_LAD.but.XL_MEF_only_peaks <- ol.ciLAD.XL.MEF_LAD$overlappingPeaks[[1]][,c(2:4)]
  overlap.with.MEF_LAD.but.XL_MEF_only_peaks <- unique(overlap.with.MEF_LAD.but.XL_MEF_only_peaks)
  output.file.name ="overlap.with.MEF_LAD.but.XL_MEF_only_peaks"
  orderPeakAndOutPut(overlap.with.MEF_LAD.but.XL_MEF_only_peaks,output.file.dir,output.file.name)
  
  overlap.with.ciLAD.but.XL_MEF_only <- toGRanges(overlap.with.ciLAD.but.XL_MEF_only_peaks)
  overlap.with.MEF_LAD.but.XL_MEF_only <- toGRanges(overlap.with.MEF_LAD.but.XL_MEF_only_peaks)
  
  ol2 <- findOverlapsOfPeaks(overlap.with.ciLAD.but.XL_MEF_only,overlap.with.MEF_LAD.but.XL_MEF_only)
  
  ol2.name <- str_sub(names(ol2$uniquePeaks),1,str_locate(names(ol2$uniquePeaks),"__")[,1]-1)
  
  overlap.with.MEF_LAD.but.XL_MEF_only.exclude.ciLAD <- ol2$uniquePeaks[ol2.name=="overlap.with.MEF_LAD.but.XL_MEF_only",]
  overlap.with.ciLAD.but.XL_MEF_only.exclude.MEF_LAD <- ol2$uniquePeaks[ol2.name=="overlap.with.ciLAD.but.XL_MEF_only",]
  
  Grange2DF <- function(overlap.with.MEF_LAD.but.XL_MEF_only.exclude.ciLAD) {
    
    gr <- overlap.with.MEF_LAD.but.XL_MEF_only.exclude.ciLAD
    
    peakAll <- data.frame(seqnames=seqnames(gr),
                          start=start(gr),
                          end=end(gr),
                          names=c(rep(".", length(gr))),
                          strands=strand(gr),peakNames=names(gr))
    peakAll 
  }
  
  overlap.with.MEF_LAD.but.XL_MEF_only.exclude.ciLAD.DF <- Grange2DF(overlap.with.MEF_LAD.but.XL_MEF_only.exclude.ciLAD)
  overlap.with.ciLAD.but.XL_MEF_only.exclude.MEF_LAD.DF <- Grange2DF(overlap.with.ciLAD.but.XL_MEF_only.exclude.MEF_LAD)
  
  output.file.name ="overlap.with.MEF_LAD.but.XL_MEF_only.exclude.ciLAD_peaks"
  orderPeakAndOutPut(overlap.with.MEF_LAD.but.XL_MEF_only.exclude.ciLAD.DF,output.file.dir,output.file.name)
  
  output.file.name ="overlap.with.ciLAD.but.XL_MEF_only.exclude.MEF_LAD_peaks"
  orderPeakAndOutPut(overlap.with.ciLAD.but.XL_MEF_only.exclude.MEF_LAD.DF,output.file.dir,output.file.name)
  
  getUniquePeaks <- function(ol.ciLAD.XL.MEF_LAD,sample.Name,output.file.dir) {
    
    ol.ciLAD.XL.MEF_LAD.uniquePeaks <- Grange2DF(ol.ciLAD.XL.MEF_LAD$uniquePeaks)
    ol.ciLAD.XL.MEF_LAD.uniquePeaks.name <- unique(str_sub(ol.ciLAD.XL.MEF_LAD.uniquePeaks$peakNames,1,str_locate(ol.ciLAD.XL.MEF_LAD.uniquePeaks$peakNames,"__")[,1]-1))
    index.sample <- which(ol.ciLAD.XL.MEF_LAD.uniquePeaks.name == sample.Name)
    XL.unique.peaks <- ol.ciLAD.XL.MEF_LAD.uniquePeaks[grep(ol.ciLAD.XL.MEF_LAD.uniquePeaks.name[index.sample],ol.ciLAD.XL.MEF_LAD.uniquePeaks$peakNames),] 
   
    output.file.name = paste0(ol.ciLAD.XL.MEF_LAD.uniquePeaks.name[index.sample],"_unique_peaks")
    output.file.name <- gsub("<","_", output.file.name)
    orderPeakAndOutPut(XL.unique.peaks,output.file.dir,output.file.name)
  }

  sample.Name <- "XL_MEF_1.5FC_minp_countF200_auto_q<e-3_chrX_q<e-2_adjPval-copy"
  getUniquePeaks(ol.ciLAD.XL.MEF_LAD,sample.Name,output.file.dir)
                             
  sample.Name <- "nonXL_MEF_1.5FC_minp_countF200_auto_q<e-3_chrX_q<e-2_adjPval-copy"
  getUniquePeaks(ol.ciLAD.non_XL.MEF_LAD,sample.Name,output.file.dir)
  
 ol.ciLAD.non_XL.MEF_LAD <- findOverlapsOfPeaks(re.out[c(1,2,3)])

  getOverLap(ol.ciLAD.non_XL.MEF_LAD,3,output.file.dir)
  getOverLap(ol.ciLAD.non_XL.MEF_LAD,2,output.file.dir)
  getOverLap(ol.ciLAD.non_XL.MEF_LAD,1,output.file.dir)
  
  overlap.with.ciLAD.but.non_XL_MEF_only_peaks <- ol.ciLAD.non_XL.MEF_LAD$overlappingPeaks[[3]][,c(11:13)]
  overlap.with.ciLAD.but.non_XL_MEF_only_peaks <- unique(overlap.with.ciLAD.but.non_XL_MEF_only_peaks)
  output.file.name ="overlap.with.ciLAD.but.non_XL_MEF_only_peaks"
  orderPeakAndOutPut(overlap.with.ciLAD.but.non_XL_MEF_only_peaks,output.file.dir,output.file.name)
  
  overlap.with.MEF_LAD.but.non_XL_MEF_only_peaks <- ol.ciLAD.non_XL.MEF_LAD$overlappingPeaks[[1]][,c(2:4)]
  overlap.with.MEF_LAD.but.non_XL_MEF_only_peaks <- unique(overlap.with.MEF_LAD.but.non_XL_MEF_only_peaks)
  output.file.name ="overlap.with.MEF_LAD.but.non_XL_MEF_only_peaks"
  orderPeakAndOutPut(overlap.with.MEF_LAD.but.non_XL_MEF_only_peaks,output.file.dir,output.file.name)
  
  overlap.with.ciLAD.but.non_XL_MEF_only <- toGRanges(overlap.with.ciLAD.but.non_XL_MEF_only_peaks)
  overlap.with.MEF_LAD.but.non_XL_MEF_only <- toGRanges(overlap.with.MEF_LAD.but.non_XL_MEF_only_peaks)
  
  ol3 <- findOverlapsOfPeaks(overlap.with.ciLAD.but.non_XL_MEF_only,overlap.with.MEF_LAD.but.non_XL_MEF_only)
  
  ol3.name <- str_sub(names(ol3$uniquePeaks),1,str_locate(names(ol3$uniquePeaks),"__")[,1]-1)
  
  overlap.with.MEF_LAD.but.non_XL_MEF_only.exclude.ciLAD <- ol3$uniquePeaks[ol3.name=="overlap.with.MEF_LAD.but.non_XL_MEF_only",]
  overlap.with.ciLAD.but.non_XL_MEF_only.exclude.MEF_LAD <- ol3$uniquePeaks[ol3.name=="overlap.with.ciLAD.but.non_XL_MEF_only",]
  
  overlap.with.MEF_LAD.but.non_XL_MEF_only.exclude.ciLAD.DF <- Grange2DF(overlap.with.MEF_LAD.but.non_XL_MEF_only.exclude.ciLAD)
  overlap.with.ciLAD.but.non_XL_MEF_only.exclude.MEF_LAD.DF <- Grange2DF(overlap.with.ciLAD.but.non_XL_MEF_only.exclude.MEF_LAD)
  
  output.file.name ="overlap.with.MEF_LAD.but.non_XL_MEF_only.exclude.ciLAD_peaks"
  orderPeakAndOutPut(overlap.with.MEF_LAD.but.non_XL_MEF_only.exclude.ciLAD.DF,output.file.dir,output.file.name)
  
  output.file.name ="overlap.with.ciLAD.but.non_XL_MEF_only.exclude.MEF_LAD_peaks"
  orderPeakAndOutPut(overlap.with.ciLAD.but.non_XL_MEF_only.exclude.MEF_LAD.DF,output.file.dir,output.file.name)
  
  outPutGrange2Bed <- function(ol.ciLAD.XL.MEF_LAD,peak.index) {
    
    gr <- ol.ciLAD.XL.MEF_LAD$peaklist[[peak.index]]
    
    peakAll <- data.frame(seqnames=seqnames(gr),
                     start=start(gr)-1,
                     end=end(gr),
                     names=c(rep(".", length(gr))),
                     strands=strand(gr),peakNames=elementMetadata(gr))
    peakAll 
  }
  
  peak1 <- outPutGrange2Bed(ol.ciLAD.XL.MEF_LAD,6)
  overlapName <- gsub("///","-overlap-",names(ol.ciLAD.XL.MEF_LAD$peaklist))[6]
  overlapName <- gsub("<","_less_",overlapName)
  orderPeakAndOutPut(peak1,output.file.dir,overlapName)
  
  peak2 <- outPutGrange2Bed(ol.ciLAD.XL.MEF_LAD,3)
  overlapName <- gsub("///","-overlap-",names(ol.ciLAD.XL.MEF_LAD$peaklist))[3]
  overlapName <- gsub("<","_less_",overlapName)
  orderPeakAndOutPut(peak2,output.file.dir,overlapName)
  
  pdf(file.path(output.file.dir,paste0("ciLAD","-XL-","MEF_LAD","_peak_overlap.pdf")))  
  venn.ciLAD.XL.MEF_LAD <- makeVennDiagram(re.out[c(1,4,3)], NameOfPeaks=c("ciLAD","XL","MEF_LAD"),
                               totalTest=5000,scaled=F,euler.d=F,by= "base",fill = c("red","blue","green"),
                               alpha = 0.50,
                               label.col = c(rep("black",7)),
                               cex = 2,
                               fontfamily = "serif",
                               fontface = "bold",
                               cat.col = c("red","blue","green"))
  
  
  venn.ciLAD.XL.MEF_LAD <- makeVennDiagram(re.out[c(1,4)], NameOfPeaks=c("ciLAD","XL"),
                                           totalTest=5000,scaled=F,euler.d=F,by= "base",fill = c("red","blue"),
                                           alpha = 0.50,
                                           label.col = c(rep("black",3)),
                                           cex = 2,
                                           fontfamily = "serif",
                                           fontface = "bold",
                                           cat.col = c("red","blue"))
  
  
  
  venn.ciLAD.XL.MEF_LAD <- makeVennDiagram(re.out[c(1,4,3)], NameOfPeaks=c("ciLAD","XL","MEF_LAD"),
                                           totalTest=10000,scaled=F,euler.d=F,by= "base",fill = c("red","blue","green"),
                                           alpha = 0.50,
                                           label.col = c(rep("black",7)),
                                           cex = 2,
                                           fontfamily = "serif",
                                           fontface = "bold",
                                           cat.col = c("red","blue","green"),connectedPeaks = "keepAll")
  
  dev.off()
  
          
  pdf(file.path(output.file.dir,paste0("XL","-MEF_LAD","_peak_overlap.pdf"))) 
  venn.XL.MEF_LAD <- makeVennDiagram(re.out[c(4,4)], NameOfPeaks=c("XL","MEF_LAD"),
                                           totalTest=10000,scaled=F, euler.d=F,fill = c("blue","green"),
                                           alpha = 0.50,
                                           label.col = c(rep("black",3)),
                                           cex = 2,
                                           fontfamily = "serif",
                                           fontface = "bold",
                                           cat.col = c("blue","green"),connectedPeaks = "keepAll")
  dev.off()
  
  u <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_9_6_2018/XL/ciLAD.mm10-bed-with-header-only-merged-in_overlap_peaks.bed"
  
  peaks=read.table(u)
  colnames(peaks)[1:3]= c("chr","start","end")
  peaks=toGRanges(peaks)

  peaks

  re.out[["ciLAD_XL"]] <- peaks
  
  ol.ciLAD_XL_MEFLAD <- findOverlapsOfPeaks(re.out[c(5,3)])
  
  overlap.with.MEF_LAD.but.non_XL_MEF_only_peaks <- ol.ciLAD_XL_MEFLAD$overlappingPeaks[[1]][,c(2:4)]
  overlap.with.MEF_LAD.but.non_XL_MEF_only_peaks <- unique(overlap.with.MEF_LAD.but.non_XL_MEF_only_peaks)
  
  peak.in.ciLAD <- ol.ciLAD_XL_MEFLAD$uniquePeaks[grep("ciLAD_XL",names(ol.ciLAD_XL_MEFLAD$uniquePeaks)),]
  
  output.file.name ="overlap.with.MEF_LAD.but.non_XL_MEF_only_peaks"
  orderPeakAndOutPut(overlap.with.MEF_LAD.but.non_XL_MEF_only_peaks,output.file.dir,output.file.name)
  
  venn.XL.MEF_LAD <- makeVennDiagram(re.out[c(5,3)], NameOfPeaks=c("MEF_LAD","ciLAD_XL"),
                                     totalTest=10000,scaled=F, euler.d=F,fill = c("blue","green"),
                                     alpha = 0.50,
                                     label.col = c(rep("black",3)),
                                     cex = 2,
                                     fontfamily = "serif",
                                     fontface = "bold",
                                     cat.col = c("blue","green"),connectedPeaks = "keepAll")
  
  ol.ciLAD_XL <- findOverlapsOfPeaks(re.out[c(1,4)],connectedPeaks = "keepAll")
  
  getCount4MergedPeaks <- function(ol.ciLAD_XL) {
    peak.name.involved.in.overlap <- ol.ciLAD_XL$mergedPeaks$peakNames
    
    X <- lapply(peak.name.involved.in.overlap, function(u){
      
      m <- length(grep("XL",u))
      n <- length(grep("ciLAD",u))
      p <- length(grep("MEF-LADs",u))
      
      x <- data.frame(XL=m,ciLAD=n,MEF_LADs=p,min_v=min(m,n,p))
      x
    })
    
    XX <- do.call(rbind.data.frame,X)
    mergedPeakId<- paste0("mergedPeak.",seq(1,dim(XX)[1]))
    
    XX <- data.frame(mergedPeakId=mergedPeakId,XX)
    XX
  }
  
  XX <- getCount4MergedPeaks(ol.ciLAD_XL)
  
  
  
  num.of.overlap.peaks.in.keepAll <- sum(XX$min_v)

  
  ol.ciLAD.XL.MEF_LAD <- findOverlapsOfPeaks(re.out[c(1,4,3)],connectedPeaks = "keepAll")
  
  YY <- getCount4MergedPeaks(ol.ciLAD.XL.MEF_LAD)
  
  ol.ciLAD.XL.MEF_LAD$mergedPeaks
  
  pdf(file.path(output.file.dir,paste0("ciLAD","-XL","_peak_overlap.pdf"))) 
  venn.ciLAD.XL <- makeVennDiagram(re.out[c(1,4)], NameOfPeaks=c("ciLAD","XL"),
                                   totalTest=10000,scaled=F,by="base", euler.d=F,fill = c("red","blue"),
                                   alpha = 0.50,
                                   label.col = c(rep("black",3)),
                                   cex = 2,
                                   fontfamily = "serif",
                                   fontface = "bold",
                                   cat.col = c("red","blue"),connectedPeaks = "keepAll")
  dev.off()
    
  
  pdf(file.path(output.file.dir,paste0("ciLAD","-non_XL-","MEF_LAD","_peak_overlap.pdf")))  
  venn.ciLAD.non_XL.MEF_LAD <- makeVennDiagram(re.out[c(1,2,3)], NameOfPeaks=c("ciLAD","non_XL","MEF_LAD"),
                                           totalTest=10000,scaled=F,euler.d=F,fill = c("red","blue","green"),
                                           alpha = 0.50,
                                           label.col = c(rep("black",7)),
                                           cex = 2,
                                           fontfamily = "serif",
                                           fontface = "bold",
                                           cat.col = c("red","blue","green"),connectedPeaks = "keepAll")
  dev.off()
  
  pdf(file.path(output.file.dir,paste0("ciLAD","-non_XL","_peak_overlap.pdf")))  
  venn.ciLAD.non_XL <- makeVennDiagram(re.out[c(1,2)], NameOfPeaks=c("ciLAD","non_XL"),
                                               totalTest=10000,scaled=F,euler.d=F,fill = c("red","blue"),
                                               alpha = 0.50,
                                               label.col = c(rep("black",3)),
                                               cex = 2,
                                               fontfamily = "serif",
                                               fontface = "bold",
                                               cat.col = c("red","blue"),connectedPeaks = "keepAll")
  dev.off()
  
  pdf(file.path(output.file.dir,paste0("non_XL-","MEF_LAD","_peak_overlap.pdf")))  
  venn.non_XL.MEF_LAD <- makeVennDiagram(re.out[c(2,3)], NameOfPeaks=c("non_XL","MEF_LAD"),
                                               totalTest=10000,scaled=F,euler.d=F,fill = c("blue","green"),
                                               alpha = 0.50,
                                               label.col = c(rep("black",3)),
                                               cex = 2,
                                               fontfamily = "serif",
                                               fontface = "bold",
                                               cat.col = c("blue","green"),connectedPeaks = "keepAll")
  
  dev.off()
  
  
  
  peak.in.ciLAD <- ol$uniquePeaks[grep("ciLAD.mm10-bed-with-header",names(ol$uniquePeaks)),]
  peak.in.MEF_LAD <- ol$uniquePeaks[-grep("ciLAD.mm10-bed-with-header",names(ol$uniquePeaks)),]
  
  pdf(file.path(output.file.dir,paste0("ciLAD-vs-","MEF_LAD","_peak_overlap.pdf")))  
  ol.by.min <- makeVennDiagram(re.out[c(1,3)], NameOfPeaks=c("ciLAD","MEF_LAD"),
                        totalTest=10000,scaled=F, euler.d=F,fill = c("red","blue"),
                        alpha = 0.50,
                        label.col = c(rep("black",3)),
                        cex = 2,
                        fontfamily = "serif",
                        fontface = "bold",
                        cat.col = c("red","blue"),connectedPeaks = "min")
  dev.off()

  re.out1 <- makeNewListAndVenn(re.out, index.pat1, index.pat2, output.file.dir, pat1, pat2)
  re.out <- makeNewListAndVenn(re.out1, index.pat1, index.pat3, output.file.dir, pat1, pat3)
  
 # genome = "Mm"
  null <- lapply(1:length(re.out),function(u,re.out,output.file.dir,genome){
    
    if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
    
    x=re.out[[u]]
    x_name=names(re.out)[u]
  
  if(genome=="Mm"){
    
    annoData <- toGRanges(EnsDb.Mmusculus.v75, feature="gene")
    
    dd.GRCm39.mm10<-toGRanges(EnsDb.Mmusculus.v75)

    overlaps.trimmed<-trim(x,use.names=TRUE)
    overlaps.anno<-annoPeaks(overlaps.trimmed,dd.GRCm39.mm10)
    write.table(overlaps.anno,file=file.path(output.file.dir,paste0(x_name,"_annotation.txt")),row.names = FALSE,quote=FALSE,sep="\t")
    
   # annotatedPeak1 <- annotatePeakInBatch(overlaps.trimmed, AnnotationData=annoData)
   #  overlaps.anno.with.entrez.id <- addGeneIDs(annotatedPeak1,"org.Mm.eg.db",IDs2Add = "symbol")
  #  write.csv(as.data.frame(unname(overlaps.anno.with.entrez.id)), file.path(output.file.dir,paste0(x_name,"_other_anno.csv")))
  #  pdf(file.path(output.file.dir,paste0(x_name,"_annotation_pie_plot.pdf")))  
  #  pie1(table(overlaps.anno.with.entrez.id$insideFeature))
  #  dev.off()
  #  write.table(annotatedPeak1,file=file.path(output.file.dir,paste0(x_name,"_annotationInBatch.txt")),row.names = FALSE,quote=FALSE,sep="\t")
    
    getGoAndPath <- function(overlaps.anno.with.entrez.id,output.file.dir,x_name) {
      over <- getEnrichedGO(overlaps.anno.with.entrez.id, orgAnn="org.Mm.eg.db",
                            maxP=0.5, minGOterm=10,
                            multiAdjMethod="BH", condense=TRUE)
      
      
      path <- getEnrichedPATH(overlaps.anno.with.entrez.id, "org.Mm.eg.db", "reactome.db", maxP=.05)
      
      write.table(path,file=file.path(output.file.dir,paste0(x_name,"_path.txt")),row.names = FALSE,quote=FALSE,sep="\t")
  
      convert2geneSymbol <- function(over1) {
        
        geneSymbol <- lapply(over1$EntrezID, function(u){
          
          x <- annotate::getSYMBOL(unlist(strsplit(as.character(u),";")), data='org.Mm.eg') 
          x
          
        })
        
        list_to_df <- function(list_for_df)
        {
          list_for_df <- as.list(list_for_df)
          
          nm <- names(list_for_df)
          if (is.null(nm)) 
            nm <- seq_along(list_for_df)
          
          df <- data.frame(name = nm, stringsAsFactors = FALSE)
          df$value <- unname(list_for_df)
          df
        }
        
        geneSymbol3 <- list_to_df(geneSymbol)
        
        over2 <- cbind(over1,geneSymbol3)
        over3 <- over2[,-which(colnames(over2) %in% c("name"))]
        colnames(over3)[which(colnames(over3) %in% c("value"))] <- "geneSymbol"
        over3
      }
      
      over$bp <- convert2geneSymbol(over$bp)
      over$mf <- convert2geneSymbol(over$mf)
      over$cc <- convert2geneSymbol(over$cc)
      
      over_bp <- as_tibble(over$bp)
      over_mf <- as_tibble(over$mf)
      over_cc <- as_tibble(over$cc)
      
      over_bp <- over_bp[order(over_bp$BH.adjusted.p.value),]
      over_mf <- over_mf[order(over_mf$BH.adjusted.p.value),]
      over_cc <- over_cc[order(over_cc$BH.adjusted.p.value),]
      
      writeTibble <- function(tibble.input, output.file.name = tempfile())
      {
        if (!dir.exists(dirname(output.file.name)))
        {
          dir.create(dirname(output.file.name), recursive = TRUE)
        }
        flatten_list = function(x)
        {
          if (typeof(x) != "list")
          {
            return(x)
          }
          sapply(x, function(y) paste(y, collapse = " ; "))
        }
        tibble.input %>% mutate_all(funs(flatten_list)) %>% write.csv(output.file.name)
      }
      
      writeTibble(over_bp, output.file.name = file.path(output.file.dir,paste0(x_name,"_GO_BP.csv")))
      writeTibble(over_mf, output.file.name = file.path(output.file.dir,paste0(x_name,"_GO_MF.csv")))
      writeTibble(over_cc, output.file.name = file.path(output.file.dir,paste0(x_name,"_GO_CC.csv")))
    }
    
    #getGoAndPath(overlaps.anno.with.entrez.id,output.file.dir,x_name)
    getGoAndPath(overlaps.anno,output.file.dir,x_name)
    
  }else if(genome=="Hs"){
    
    library(EnsDb.Hsapiens.v75)
    #annoData<-toGRanges(EnsDb.Hsapiens.v75, feature="gene")
    
    dd.hs<-toGRanges(EnsDb.Hsapiens.v75)
    
    print(seqinfo(dd.hs))
    print(seqlevels(dd.hs))
    
    #print(seqlevels(dd.hs)[,1])
    
    #print(seqlevels(re.out[[1]])[,1])
    
    # seqlevels(dd.hs,force=TRUE) <- c("chr1","chr10","chr11","chr12","chr13",
    #                                           "chr14","chr15","chr16","chr17","chr18","chr19","chr2",
    #                                           "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY")
    
    #temp4=
    
    re.out.L<-lapply(1:length(re.out),function(u,re.out,dd.hs){
      
      x=re.out[[u]]
      x_name=names(re.out)[u]
      
      seqlevels(dd.hs,force=TRUE)<-seqinfo(x)@seqnames
      #print(seqinfo(re.out.trimmed))
      #print(seqlevels(re.out.trimmed))
      seqinfo(x)<-seqinfo(dd.hs)
      #GRCm38/mm10
      #dd<-toGRanges(EnsDb.Mmusculus.v79)
      #seqinfo(dd)
      #library(ensembldb)
      #library(GenomeInfoDb)
      seqlevelsStyle(x) <- seqlevelsStyle(dd.hs)
      re.out.trimmed<-trim(x, use.names=TRUE)
      overlaps.anno<-annoPeaks(re.out.trimmed,dd.hs)
      
      write.table(overlaps.anno,file=paste0(temp3,"/",x_name,"_annotation.txt"),row.names = FALSE,quote=FALSE,sep="\t")
    },re.out,dd.hs)
    
  }
  
  },re.out,output.file.dir,genome)
  
  
  peaks1 <- GRanges(seqnames=c("1", "1", "1"),
                    IRanges(start=c(8, 15, 60),
                            end=c(20, 50, 100), 
                            names=c("Site1", "Site2", "Site3")),
                            strand="+",
                    feature=c("a","a","a"))
  
  peaks2 = GRanges(seqnames=c("1", "2", "3", "1", "2"), 
                   IRanges(start = c(10, 30,90,
                                     100,200),
                           end = c(20,50, 200, 
                                   150, 500),
                           names = c("t1", "t2", "t3", "t4", "t5")), 
                   strand = c("+", "+", "+", "+", "+"), 
                   feature=c("a","a","a","a","a"))
  
   venn.plot <- venn.diagram(
         x = list(
               A = c(1:10),
               B = c(11:90),
               C = c(81:90)
           ),
         euler.d = TRUE,
         scaled = TRUE, 
         filename = "Euler_3set_scaled.tiff",
         cex = 2.5,
         cat.cex = 2.5,
         cat.pos = 0);
  
  peaks1 <- GRanges(seqnames=c("1"),
                    IRanges(start=c(2),
                            end=c(10), 
                            names=c("Site1")),
                    strand="+",
                    feature=c("a"))
  peaks2 = GRanges(seqnames=c("1"), 
                   IRanges(start = c(8,30),
                           end = c(20,40),
                           names = c("t1","t2")), 
                   strand = c("+","+"), 
                   feature=c("a","a"))
  
  peaks1 <- re.out[[1]]
  peaks2 <- re.out[[4]]
  
  getCount4Venn <- function(re.out,peak.index,name,output.file.dir) {
    
    name <- c("F121-9","XL_MEFs")
    
    peak.index <- c(1,2)
    grl <- GRangesList(re.out[peak.index])
    names(grl) <- name
    
    Z <- ChIPpeakAnno:::vennCounts(grl,length(peak.index),names=names(grl),by="base") 
    
    if(length(peak.index)==2){
    ZZ <- Z$venn_cnt[-which(row.names(Z$venn_cnt)=="00"),]
    }
    
    if(length(peak.index)==3){
      ZZ <- Z$venn_cnt[-which(row.names(Z$venn_cnt)=="000"),]
    }
    
    y <- ZZ[,"Counts"]
    names(y) <- row.names(ZZ)
    
    labels = name
    pdf(file = file.path(output.file.dir,paste0(paste(name,collapse = "-"),".pdf")))
    plot.new()
    colorfulVennPlot::plotVenn(y, labels, Colors=rainbow(7))
    dev.off()
    
    fit1 <- euler(c("ciLAD" = 785, "LAD" = 565, "nonXL_MEF" = 167,
                    "ciLAD&LAD" = 3, "ciLAD&nonXL_MEF" = 101, "LAD&nonXL_MEF" = 541,
                    "ciLAD&LAD&nonXL_MEF" = 2),shape = "ellipse")
    
    ZZZZ <- round(ZZZ/1000000,0)
    
    fit2 <- euler(ZZZZ,shape = "ellipse")
    fit2.plot <- plot(fit2,quantities = TRUE,fill = rainbow(7),lty = 1:2,labels = list(font = 1),alpha=0.7)
    title <- textGrob("F121_9-vs-XL_MEF", x= unit(250,"native"),y= unit(30,"native"))
    footnote <- textGrob("Numbers are nucleotides in MB", x=unit(380,"native"),y=unit(500,"native"))
    gt <- gTree(children=gList(fit2.plot, title, footnote))
    gt2 <- grid.draw(gt)
    
    grid.newpage()
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
    
    png("g.png",width = 1000, height = 600, units = "px") 
    grid.draw(g) 
    dev.off()
    
    png("F121_9-vs-XL_MEF_overlapping.png",width = 1000, height = 600)
    grid.newpage()
    grid.draw(gt)
    dev.off()
    
    grid.arrange(gt,ncol=1,nrow=1)
    ggsave("saving.pdf")

    pdf(file="filename.pdf", onefile=FALSE) # or other device
   
    dev.off()
    
    grid.ls(gt)
    
    t <- grid.get("GRID.text.73")
    t$label <- "F121-9-vs-XLMEF"
    t$x <- unit(250,"native")
    t$y <- unit(40,"native")
    grid.edit("GRID.text.73",label=t$label, x=t$x,y=t$y)
    
    s <- grid.get("GRID.text.74")
    s$label <- "Numbers are nucleotides in MB"
    s$x <- unit(250,"native")
    s$y <- unit(10,"native")
    grid.edit("GRID.text.74",label=s$label,x=s$x,s=s$y)
    
    
    
    
    colorfulVennPlot::plotVenn(ZZ, labels=names(ZZZZ)[1:3], Colors=rainbow(7))
    
    fit1 <- euler(c("ciLAD" = 785, "LAD" = 565, "nonXL_MEF" = 167,
                    "ciLAD&LAD" = 3, "ciLAD&nonXL_MEF" = 101, "LAD&nonXL_MEF" = 541,
                    "ciLAD&LAD&nonXL_MEF" = 2),shape = "ellipse")
    
    grid.newpage()
    fit1.plot <- plot(fit1,quantities = TRUE,fill = rainbow(7),lty = 1:2,labels = list(font = 1),alpha=0.7)
    fit1.plot
    grid.ls(fit1.plot)
    t <- grid.get("quantities.grob")
    names(t)
    t$label
    t$x
    t$y
    
    ChangeXY4CertainLabel <- function(x,which.label){
      options(digits=15)
      Y <- lapply(1:length(x),function(u,x,which.label){
        
        if(u == which.label){
          cat("old_value is",x[u],"\n")
          
          cat('new value you want:', '\n')
          
          m<-scan("",n=1, quiet=TRUE)
          
          z = m 
        }else
        {
          z=x[u]
        }
        
        z      
        
      },x,which.label)
      
      Y <- unit(as.array(unlist(Y)),"native")
      Y
    }
    
    new.y <- ChangeXY4CertainLabel(t$y,4)
    grid.edit("quantities.grob",y=new.y)
    
    dev.off()
  
    y
  }
  
  peak.index <-c(1,4)
  name <- c("ciLAD","XL_MEF") 
  YY <- getCount4Venn(re.out,peak.index,name,output.file.dir) 
  
  peak.index <-c(1,4,3)
  name <- c("ciLAD","XL_MEF","LAD") 
  

  g1 <- unique(combn(name,1,simplify=FALSE))
  g2 <- unique(combn(name,2,simplify=FALSE))
  g3 <- unique(combn(name,3,simplify=FALSE))
  
  g2 <- unlist(lapply(g2, function(u){paste(u,collapse = "&")}))
  g3 <- unlist(lapply(g3, function(u){paste(u,collapse = "&")}))
  
  LabelName <- unlist(c(g1,g2,g3))
  
  LabelName <- unlist(c(g1,g2))
  
  xx <- c(4,2,1,6,5,3,7)
  xx <- c(2,1,3)
  names(xx) <- LabelName
  
  ZZ <- unlist(lapply(xx, function(u){
    z <- paste0(as.binary(u,n=2), collapse = "")
    z
  }))
  
  ZZZ <- Z$venn_cnt[match(ZZ,row.names(Z$venn_cnt)),'Counts']
  names(ZZZ) <- names(ZZ)
  
  YY <- getCount4Venn(re.out,peak.index,name,output.file.dir) 
  
  peak.index <-c(1,2,3)
  name <- c("ciLAD","nonXL_MEF","LAD") 
  YY <- getCount4Venn(re.out,peak.index,name,output.file.dir) 
  
  labels <- names(grl)
  
  plot.new()
  colorfulVennPlot::plotVenn(y, labels, Colors=rainbow(7))
  
  
  test.ol.2 <- makeVennDiagram(list(peaks1, peaks2), NameOfPeaks=c("TF1", "TF2"),
                               totalTest=100000,scaled=TRUE, euler.d=TRUE, by="base")
 
  test.ol.2 <- makeVennDiagram(list(re.out[[1]],re.out[[4]],re.out[[3]]), NameOfPeaks=c("ciLAD", "XL_MEF","MEF_LAD"),
                               totalTest=100000,scaled=TRUE, euler.d=TRUE, by="base")
  
  test.ol.2 <- makeVennDiagram(list(re.out[[1]],re.out[[4]]), NameOfPeaks=c("ciLAD", "XL_MEF"),
                               totalTest=100000,scaled=TRUE, euler.d=TRUE, by="base")
  
  test.ol.2$vennCounts
  
}

#' AnnotatePeak2
#'
#' @param input.file.dir 
#' @param input.file.pattern 
#' @param index.file 
#' @param output.file.dir 
#' @param genome 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' input.file.dir="/media/H_driver/2016/Danny/Danny_chip_PeakCall/"
#' input.file.pattern="*macs142_peaks.bed"
#' output.file.dir="/media/H_driver/2016/Danny/Danny_chip_PeakCall/"
#'
#' AnnotatePeak(input.file.dir,input.file.pattern,8,output.file.dir,genome="Mm")
#' AnnotatePeak(input.file.dir,input.file.pattern,7,output.file.dir,genome="Hs")
#' 
#' 
#' 
AnnotatePeak2 <- function(input.file.dir,input.file.pattern,index.file,output.file.dir,genome) {
  
  library(ChIPpeakAnno)
  
  dir.name=input.file.dir
  input.file.pattern=input.file.pattern
  
  dir.name=reformatPath(dir.name)
  output.dir.name=reformatPath(output.file.dir)
  
  #print(output.dir.name)
  #temp=Sys.time()
  #temp1=gsub(":","-",Sys.time())
  #temp2=gsub(" ","-",temp1)
  #temp3=paste0(output.dir.name,"AnalysisResults_at_",temp2)
  
  temp3=output.dir.name
  
  dir.create(temp3)
  
  file.name=paste0(dir.name,dir(dir.name,recursive = TRUE,pattern=input.file.pattern))
  file.name.2<-as.list(file.name)
  
  names(file.name.2)=sapply(strsplit(file.name,split="\\/"),"[[",index.file)
  
  print(file.name.2)
  
  file.name.3<-file.name.2
  
  #sample.name<-sapply(strsplit(names(file.name.3),split="_peaks_"),"[[",1)
  
  #names(file.name.3)=sample.name
  
  file.name.4 <-file.name.3
  
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
  
  if(genome=="Mm"){
    
    annoData <- toGRanges(EnsDb.Mmusculus.v75, feature="gene")
    
    ol <- findOverlapsOfPeaks(re.out[c(2,4,1)])
    
    overlaps<-ol$peaklist$`11_2470IUPUI_WT_BM_SMC1_peaks.bed///13_2470IUPUI_WT_BM_Rad21_peaks.bed///10_WT_BM_ASXL1_peaks.bed`
    
    binOverFeature(overlaps, annotationData=annoData,
                   radius=5000, nbins=20, FUN=length, errFun=0,
                   ylab="count",
                   main="Distribution of aggregated peak numbers around TSS")
    
    overlaps.trimmed<-trim(overlaps, use.names=TRUE)
    
    library(EnsDb.Mmusculus.v79)
    dd.GRCm39.mm10<-toGRanges(EnsDb.Mmusculus.v75)
    #seqinfo(dd.GRCm39.mm10)
    #seqlevels(dd.GRCm39.mm10)
    
    seqlevels(dd.GRCm39.mm10,force=TRUE) <- c("chr1","chr10","chr11","chr12","chr13",
                                              "chr14","chr15","chr16","chr17","chr18","chr19","chr2",
                                              "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY")
    
    seqinfo(overlaps)<-seqinfo(dd.GRCm39.mm10)
    #GRCm38/mm10
    #dd<-toGRanges(EnsDb.Mmusculus.v79)
    #seqinfo(dd)
    #library(ensembldb)
    #library(GenomeInfoDb)
    seqlevelsStyle(overlaps.trimmed) <- seqlevelsStyle(dd.GRCm39.mm10)
    overlaps.anno<-annoPeaks(overlaps.trimmed,dd.GRCm39.mm10)
    write.table(overlaps.anno,file=paste0(temp3,"/","annotation.txt"),row.names = FALSE,quote=FALSE,sep="\t")
  }else if(genome=="Hs"){
    
    library(EnsDb.Hsapiens.v75)
    #annoData<-toGRanges(EnsDb.Hsapiens.v75, feature="gene")
    
    dd.hs<-toGRanges(EnsDb.Hsapiens.v75)
    
    print(seqinfo(dd.hs))
    print(seqlevels(dd.hs))
    
    #print(seqlevels(dd.hs)[,1])
    
    #print(seqlevels(re.out[[1]])[,1])
    
    # seqlevels(dd.hs,force=TRUE) <- c("chr1","chr10","chr11","chr12","chr13",
    #                                           "chr14","chr15","chr16","chr17","chr18","chr19","chr2",
    #                                           "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY")
    
    #temp4=
    
    re.out.L<-lapply(1:length(re.out),function(u,re.out,dd.hs){
      
      x=re.out[[u]]
      x_name=names(re.out)[u]
      
      seqlevels(dd.hs,force=TRUE)<-seqinfo(x)@seqnames
      #print(seqinfo(re.out.trimmed))
      #print(seqlevels(re.out.trimmed))
      seqinfo(x)<-seqinfo(dd.hs)
      #GRCm38/mm10
      #dd<-toGRanges(EnsDb.Mmusculus.v79)
      #seqinfo(dd)
      #library(ensembldb)
      #library(GenomeInfoDb)
      seqlevelsStyle(x) <- seqlevelsStyle(dd.hs)
      re.out.trimmed<-trim(x, use.names=TRUE)
      overlaps.anno<-annoPeaks(re.out.trimmed,dd.hs)
      
      write.table(overlaps.anno,file=paste0(temp3,"/",x_name,"_annotation.txt"),row.names = FALSE,quote=FALSE,sep="\t")
    },re.out,dd.hs)
    
  }
  
}

#' AnnotatePeak3
#'
#' @param input.file.dir 
#' @param output.file.dir 
#' @param genome 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' input.file.dir="/media/H_driver/2016/Danny/Danny_chip_PeakCall/"
#' output.file.dir="/media/H_driver/2016/Danny/Danny_chip_PeakCall/"
#'
#' AnnotatePeak3(input.file.dir,output.file.dir,genome="Mm")
#' AnnotatePeak3(input.file.dir,output.file.dir,genome="Hs")
#' 
#' 
#' 
AnnotatePeak3<- function(input.file.dir,output.file.dir,genome) {
  
  library(ChIPpeakAnno)
  
  re<-ParserReadFiles(input.file.dir,"bed",output.file.dir)
  
  #print(file.name.2)
  
  file.name.2<-re$input
  output.dir.name=output.file.dir
  
  #dir.name=input.file.dir
  #input.file.pattern=input.file.pattern
  
  #dir.name=reformatPath(dir.name)
  
  #output.dir.name=reformatPath(output.file.dir)
  
  #print(output.dir.name)
  #temp=Sys.time()
  #temp1=gsub(":","-",Sys.time())
  #temp2=gsub(" ","-",temp1)
  #temp3=paste0(output.dir.name,"AnalysisResults_at_",temp2)
  
  temp3=output.dir.name
  
  dir.create(temp3)
  
  #file.name=paste0(dir.name,dir(dir.name,recursive = TRUE,pattern=input.file.pattern))
  #file.name.2<-as.list(file.name)
  
  #names(file.name.2)=sapply(strsplit(file.name,split="\\/"),"[[",index.file)
  
  #print(file.name.2)
  
  file.name.3<-file.name.2
  
  #sample.name<-sapply(strsplit(names(file.name.3),split="_peaks_"),"[[",1)
  
  #names(file.name.3)=sample.name
  
  file.name.4 <-file.name.3
  
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
  
  if(genome=="Mm"){
    
    annoData <- toGRanges(EnsDb.Mmusculus.v75, feature="gene")
    
    ol <- findOverlapsOfPeaks(re.out[c(2,4,1)])
    
    overlaps<-ol$peaklist$`11_2470IUPUI_WT_BM_SMC1_peaks.bed///13_2470IUPUI_WT_BM_Rad21_peaks.bed///10_WT_BM_ASXL1_peaks.bed`
    
    binOverFeature(overlaps, annotationData=annoData,
                   radius=5000, nbins=20, FUN=length, errFun=0,
                   ylab="count",
                   main="Distribution of aggregated peak numbers around TSS")
    
    overlaps.trimmed<-trim(overlaps, use.names=TRUE)
    
    library(EnsDb.Mmusculus.v79)
    dd.GRCm39.mm10<-toGRanges(EnsDb.Mmusculus.v75)
    #seqinfo(dd.GRCm39.mm10)
    #seqlevels(dd.GRCm39.mm10)
    
    seqlevels(dd.GRCm39.mm10,force=TRUE) <- c("chr1","chr10","chr11","chr12","chr13",
                                              "chr14","chr15","chr16","chr17","chr18","chr19","chr2",
                                              "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY")
    
    seqinfo(overlaps)<-seqinfo(dd.GRCm39.mm10)
    #GRCm38/mm10
    #dd<-toGRanges(EnsDb.Mmusculus.v79)
    #seqinfo(dd)
    #library(ensembldb)
    #library(GenomeInfoDb)
    seqlevelsStyle(overlaps.trimmed) <- seqlevelsStyle(dd.GRCm39.mm10)
    overlaps.anno<-annoPeaks(overlaps.trimmed,dd.GRCm39.mm10)
    write.table(overlaps.anno,file=paste0(temp3,"/","annotation.txt"),row.names = FALSE,quote=FALSE,sep="\t")
  }else if(genome=="Hs"){
    
    library(EnsDb.Hsapiens.v75)
    #annoData<-toGRanges(EnsDb.Hsapiens.v75, feature="gene")
    
    dd.hs<-toGRanges(EnsDb.Hsapiens.v75)
    
    print(seqinfo(dd.hs))
    print(seqlevels(dd.hs))
    
    #print(seqlevels(dd.hs)[,1])
    
    #print(seqlevels(re.out[[1]])[,1])
    
    # seqlevels(dd.hs,force=TRUE) <- c("chr1","chr10","chr11","chr12","chr13",
    #                                           "chr14","chr15","chr16","chr17","chr18","chr19","chr2",
    #                                           "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY")
    
    #temp4=
    
    re.out.L<-lapply(1:length(re.out),function(u,re.out,dd.hs){
      
      x=re.out[[u]]
      x_name=names(re.out)[u]
      
      seqlevels(dd.hs,force=TRUE)<-seqinfo(x)@seqnames
      #print(seqinfo(re.out.trimmed))
      #print(seqlevels(re.out.trimmed))
      seqinfo(x)<-seqinfo(dd.hs)
      #GRCm38/mm10
      #dd<-toGRanges(EnsDb.Mmusculus.v79)
      #seqinfo(dd)
      #library(ensembldb)
      #library(GenomeInfoDb)
      seqlevelsStyle(x) <- seqlevelsStyle(dd.hs)
      re.out.trimmed<-trim(x, use.names=TRUE)
      overlaps.anno<-annoPeaks(re.out.trimmed,dd.hs)
      
      write.table(overlaps.anno,file=paste0(temp3,"/",x_name,"_annotation.txt"),row.names = FALSE,quote=FALSE,sep="\t")
    },re.out,dd.hs)
    
  }
  
}

#' AnntationUsingChipSeeker
#' 
#' Use ChipSeeker to annotate data
#'
#' @param dir.name: path for bed files
#' @param input.file.pattern : input file pattern
#' @param out.dir.name output file directory
#' @param txdb Annotation databased to be used("hg19","hg38")
#' @param DD distance definition around TSS
#' 
#' @return
#' @export
#'
#' @examples
#' 
#' dir.name="/Volumes/Bioinformatics$/2017/DannyNewData/BindDiff/common_peaks_bed"
#' input.file.pattern="common_peaks.bed"
#' out.dir.name="/Volumes/Bioinformatics$/2017/DannyNewData/AnnotationNew4DBA"
#' txdb="hg19"
#' DD=5000
#' 
#' AnntationUsingChipSeeker(dir.name,input.file.pattern,out.dir.name,txdb=txdb,DD,distanceToTSS_cutoff=5000, AP=c("Promoter","Intron"))
#'
#' res.promoter <- AnntationUsingChipSeeker(dir.name,input.file.pattern,out.dir.name,txdb=txdb,DD,distanceToTSS_cutoff=5000,AP=c("Promoter"))
#' 
#' AnntationUsingChipSeeker(dir.name,input.file.pattern,out.dir.name,txdb=txdb,DD,distanceToTSS_cutoff=5000,AP=c("Intron"))
#'   
AnntationUsingChipSeeker <- function(dir.name,input.file.pattern,out.dir.name,txdb=c("hg19","hg38"),DD,distanceToTSS_cutoff=5000,assignGenomicAnnotation=TRUE,AP=c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic")) {
  
  #re<-ParserReadFiles(dir.name,input.file.pattern)
  
  #re.bed<-re$input
  
  #re.peaks.only.bed.2 <- re.bed
  
  file.1 <- list.files(dir.name,pattern=input.file.pattern,all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
  
  #re.bed<-re$input
  
  re.peaks.only.bed.2 <- file.1
  
  # if(length(dir(dir.name,pattern="peaks.bed"))!=0)
  # {
  # re.peaks.only.bed.2<-FL(re.bed,'peaks')
  # cat("peaks\n")
  # print(re.peaks.only.bed.2)
  # }
  # 
  # if(length(dir(dir.name,pattern="summits.bed"))!=0){
  # re.summits.only.bed<-FL(re.bed,'summits')
  # cat("summits\n")
  # print(re.summits.only.bed)
  # }
  
  txdb<-match.arg(txdb)
  
  switch (txdb,
          hg38 = {
            cat("Use hg38\n")
            txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
          },
          {
            cat("Use hg19\n") 
            txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
          }
  )
  
  APpath <- paste(AP,collapse = "_")
  
  temp3=file.path(out.dir.name,"Annotation",APpath)
  
  if(!dir.exists(temp3)){dir.create(temp3,recursive = TRUE)}
  
  d=DD
  res <- lapply(1:length(re.peaks.only.bed.2),function(u,re.peaks.only.bed.2,d){
    
    peaks=readPeakFile(re.peaks.only.bed.2[[u]],as="data.frame")
    
    print(head(peaks))
    
    peakAnno <- annotatePeak(re.peaks.only.bed.2[[u]], tssRegion=c(-d, d),
                             TxDb=txdb,assignGenomicAnnotation=assignGenomicAnnotation,genomicAnnotationPriority=AP,annoDb="org.Hs.eg.db")
    
    #select.index <- which(peakAnno$distanceToTSS<20,000 && peakAnno$distanceToTSS >= -20,000)
    dropAnnoM <- function (csAnno, distanceToTSS_cutoff) 
    {
      idx <- which(abs(mcols(csAnno@anno)[["distanceToTSS"]]) < 
                     distanceToTSS_cutoff)
      csAnno@anno <- csAnno@anno[idx]
      csAnno@peakNum <- length(idx)
      if (csAnno@hasGenomicAnnotation) {
        csAnno@annoStat <- ChIPseeker:::getGenomicAnnoStat(csAnno@anno)
        csAnno@detailGenomicAnnotation = csAnno@detailGenomicAnnotation[idx, 
                                                                        ]
      }
      csAnno
    }
    
    peakAnno <- dropAnnoM(peakAnno,distanceToTSS_cutoff = distanceToTSS_cutoff)
    
    x_name=names(re.peaks.only.bed.2)[u]
    cat(x_name)
    png(file.path(temp3,paste0(x_name,"_",d,APpath,"_around_tss_annotation_pie.png")))
    plotAnnoPie(peakAnno)
    dev.off()
    
    peaks.anno=as.data.frame(peakAnno)
    
    print(head(peaks.anno))
    
    #print(paste0(peaks[,c(2,3)]))
    print(colnames(peaks.anno))
    write.table(peaks.anno,file=file.path(temp3,paste0(x_name,"_",d,APpath,"_around_tss_annotation_4_only_mapped_peaks.xls")),
                row.names = FALSE,quote=FALSE,sep="\t")
    
    unmapped.peaks<-peaks[-which(paste0(peaks[,2],"_",peaks[,3]) %in% paste0(peaks.anno[,2],"_",peaks.anno[,3])),]
    
    cat(dim(peaks)[1]," ",dim(peaks.anno)[1]," ",dim(unmapped.peaks)[1],"\n")
    
    
    if(dim(unmapped.peaks)[1]!=0){
      
      colnames(unmapped.peaks)=colnames(peaks.anno)[1:6]
      
      unmapped.peaks.3<-smartbind(peaks.anno,unmapped.peaks)
      
      unmapped.peaks.4<-unmapped.peaks.3[order(unmapped.peaks.3[,1],unmapped.peaks.3[,2]),]
      
      write.table(unmapped.peaks.4,file=file.path(temp3,paste0(x_name,"_",d,APpath,"_around_tss_annotation_4_all_peaks.xls")),row.names = FALSE,quote=FALSE,sep="\t")
    }
    
    re <- list(peaksMappedOnly=peaks.anno,peaksAll=unmapped.peaks.4)
    re
  },re.peaks.only.bed.2,d)
  
  res
}

#' BamFileSortIndexVisualization 
#'
#' @param input.file.dir 
#' @param output.file.dir 
#' @param genome 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' input.file.dir="/projects/scratch/bbc/Project/Danny_chip/Alignment/BWA/"
#' output.file.dir="/scratch/projects/bbc/aiminy_project/"
#' genome="Hs"
#' d=4000
#' BamFileSortIndexVisualization(input.file.dir,output.file.dir,4000,genome)
#' 
#'bsub -P bbc -J "zhaoCJun" -o %J.zhaoCJun.log -e %J.zhaoCJun.err -W 72:00 -n 32 -q parallel -R 'rusage[mem= 16000 ] span[ptile= 16 ]' -u aimin.yan@med.miami.edu R -e 'library(ChipSeq);re <- ChipSeq:::BamFileSortIndexVisualization("/projects/scratch/bbc/Project/Danny_chip/Alignment/BWA","/scratch/projects/bbc/aiminy_project/7_27_2017",400000,"Hs")'
#'
#' 
BamFileSortIndexVisualization <- function(input.file.dir,input.pattern,output.file.dir,d,genome,bam.sort=FALSE) {
  
  #library(ChIPpeakAnno)
  
  #re<-ParserReadFiles(input.file.dir,"bam")

  file.1 <- list.files(input.file.dir,pattern=input.pattern, all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
  
  file.name.2<-file.1
  #output.dir.name=re$output
  
  #temp3=paste0(output.dir.name,"_visualization")
  
  temp3=output.file.dir
  
  if (!dir.exists(temp3))
  {
    dir.create(temp3, recursive = TRUE)
  }
    
  #dir.create(temp3)
  
  re.out<-file.name.2
  
  if(bam.sort==FALSE){
  
  cmd1="samtools sort"
  
  lapply(1:length(re.out),function(u,re.out,temp3){
    
    x=re.out[[u]]
    x_name=tools::file_path_sans_ext(basename(x))
    cmd2=paste0(cmd1," ",x," ","-o"," ",file.path(temp3,paste0(x_name,"_sorted.bam")))
    
    print(cmd2)
    system(cmd2)
    
  },re.out,temp3)
  
  cmd3="samtools index"
  
  lapply(1:length(re.out),function(u,re.out,temp3){
    
    x=re.out[[u]]
    x_name=tools::file_path_sans_ext(basename(x))
    cmd4=paste0(cmd3," ",file.path(temp3,paste0(x_name,"_sorted.bam")))
    
    #print(cmd2)
    system(cmd4)
    
  },re.out,temp3)
}
  
  cmd5="ngs.plot.r -G hg19 -R tss -C"
  cmd6="-O"
  #cmd7="-L 4000"
  cmd7=paste("-L",d,sep=" ")
  
  #cmd3="-L 4000 -RR 1 -CD 1 -CO \\\"blue\\\""
  
  #ngs.plot.r -G hg19 -R tss -C $1 -O $2 -L 4000 -RR 1 -CD 1 -CO "blue"
  
  #file.name.3<-file.name.2[-6]
  
  lapply(1:length(re.out),function(u,re.out,temp3){
    
    x=re.out[[u]]
    x_name=tools::file_path_sans_ext(basename(x))
    
    cmd8=paste(cmd5,file.path(temp3,paste0(x_name,"_sorted.bam")),cmd6,file.path(temp3,x_name),cmd7,sep=" ")
    
    print(cmd8)
    system(cmd8, intern = TRUE, ignore.stderr = TRUE)
    
    #re=read.table(u,header=FALSE)
    #  re<-as.character(re[,1])
    #  #colnames(re)=c("Count","GeneName")
    #  re
  },re.out,temp3)
  
  
}

convertBam2StrandBw2 <- function(input.bam.file.dir, output.bw.file.dir, BigMem = FALSE,
                                 cores = 15, Memory = 25000, Wall.time = "72:00", span.ptile = 8)
{
  re <- parserreadfiles(input.bam.file.dir, "bam")
  
  res <- re$input
  
  m.id <- grep("login", system("hostname", intern = TRUE))
  
  if (!dir.exists(output.bw.file.dir))
  {
    dir.create(output.bw.file.dir, recursive = TRUE)
  }
  
  # job.name=paste0('bamSort[',length(res),']')
  
  cmd.l <- lapply(1:length(res), function(u, m.id, Wall.time, cores, Memory,
                                          span.ptile, res, output.bw.file.dir)
  {
    # path_name = dirname(res[[u]]) path_name2 <- basename(path_name)
    
    file_name = file_path_sans_ext(basename(res[[u]]))
    
    # file_name <- paste0(path_name2,'-',file_name)
    #u <- 3
    if (m.id == 1)
    {
      if (BigMem == TRUE)
      {
        cmd0 = paste(Wall.time, "-n", cores, "-q bigmem -R 'rusage[mem=",
                     Memory, "] span[ptile=", span.ptile, "]' -u aimin.yan@med.miami.edu",
                     sep = " ")
      } else
      {
        cmd0 = paste(Wall.time, "-n", cores, "-q general -R 'rusage[mem=",
                     Memory, "] span[ptile=", span.ptile, "]' aimin.yan@med.miami.edu",
                     sep = " ")
      }
      
      job.name = paste0("bam2wig.", u)
      cmd1 = paste0("bsub -P bbc -J \"", job.name, paste0("\" -o %J.",
                                                          job.name, ".log "), paste0("-e %J.", job.name, ".err -W"))
      
      # job.name=paste0('bamSort[',length(res),']') cmd1 = paste0('bsub -w
      # \'done(\'bamSort[*]\')\'', 'bsub -P bbc -J \'',job.name,paste0('\'
      # -o %J.log '),paste0('-e %J.err -W')) job.name=paste0('Bdg[',u,']') cmd1 =
      # paste0('bsub -w \'done(\'bamIndex[*]\') && done(\'Chrosome\')\'',
      # 'bsub -P bbc -J \'',job.name,paste0('\' -o %J.',job.name,'.log
      # '),paste0('-e %J.',job.name,'.err -W'))
      if (u <= 6)
      {
        cmd2 = paste("bam2wig.pl -pe --pos span --strand --bw  --bwapp $HOME/kentUtils/bin/wigToBigWig --out",
                     file.path(output.bw.file.dir, paste0(file_name, "_2.bw")),
                     "--in", res[[u]], sep = " ")
      } else
      {
        cmd2 = paste("bam2wig.pl --pos span --strand --bw --bwapp $HOME/kentUtils/bin/wigToBigWig --out",
                     file.path(output.bw.file.dir, paste0(file_name, "_2.bw")),
                     "--in", res[[u]], sep = " ")
      }
      cmd3 = paste(cmd1, cmd0, cmd2, sep = " ")
    } else
    {
      if (u <= 6)
      {
        cmd3 = paste("bam2wig.pl -pe --pos span --strand --bw --bwapp $HOME/kentUtils/bin/wigToBigWig --out",
                     file.path(output.bw.file.dir, paste0(file_name, "_2.bw")),
                     "--in", res[[u]], sep = " ")
      } else
      {
        cmd3 = paste("bam2wig.pl --pos span --strand --bw --bwapp $HOME/kentUtils/bin/wigToBigWig --out",
                     file.path(output.bw.file.dir, paste0(file_name, "_2.bw")),
                     "--in", res[[u]], sep = " ")
      }
    }
    
    cmd <- cmd3
    
    system(cmd)
    
    cmd
  }, m.id, Wall.time, cores, Memory, span.ptile, res, output.bw.file.dir)
  
}

#'R -e 'library(ChipSeq);ChipSeq:::plotBam(input.file.dir="/scratch/projects/bbc/Project/Danny_chip2/Alignment/BWA",file.type="*marked.bam",output.file.dir="/scratch/projects/bbc/aiminy_project/DannyNewNgsPlot")'

#'R -e 'library(ChipSeq); x <- ChipSeq:::plotBam(input.file.dir="/scratch/projects/bbc/Project/Danny_chip2/Alignment/BWA",file.type="*marked.bam",output.file.dir="/scratch/projects/bbc/aiminy_project/DannyNewNgsPlot",cores = 8, Memory = 16000,span.ptile = 4,wait = FALSE)'

#'R -e 'library(ChipSeq); x <- ChipSeq:::plotBam(input.file.dir="/scratch/projects/bbc/Project/Danny_chip2/Alignment/BWA",file.type="*marked.bam",output.file.dir="/scratch/projects/bbc/aiminy_project/DannyNewNgsPlot",job.option = "parallel", cores = 32, Memory = 16000,span.ptile = 16,wait = TRUE)'

plotBam <- function(input.file.dir,file.type,output.file.dir,job.option = "general",cores = 15, Memory = 25000, Wall.time = "72:00", span.ptile = 8,wait=TRUE) {
  
  #library(ChIPpeakAnno)
  
  #cat(job.option,"\n")
  #job.option <- as.character(job.option)
  #cat(job.option,"\n")
  
  re<-ParserReadFiles(input.file.dir,file.type)
  
  file.name.2<-re$input
  #output.dir.name=re$output
  
  re.out<-file.name.2
  
  m.id <- grep("login", system("hostname", intern = TRUE))
  
  if (!dir.exists(output.file.dir))
  {
    dir.create(output.file.dir, recursive = TRUE)
  }
  
  temp3= output.file.dir
  
  #cmd1="samtools sort"
  
  cmd.1 <- lapply(1:length(re.out),function(u,m.id,job.option,Wall.time, cores, Memory,
                                            span.ptile,re.out,temp3){
    
    file_name = file_path_sans_ext(basename(re.out[[u]]))
    
    # file_name <- paste0(path_name2,'-',file_name)
    #u <- 3
    if (m.id == 1)
    {
      
      #job.option <- match.arg(job.option)
      
      switch (job.option,
              parallel = {
                cmd0 = paste(Wall.time, "-n", cores, "-q parallel -R 'rusage[mem=",
                             Memory, "] span[ptile=", span.ptile, "]' -u aimin.yan@med.miami.edu",
                             sep = " ")        
              },
              bigmem = {
                cmd0 = paste(Wall.time, "-n", cores, "-q bigmem -R 'rusage[mem=",
                             Memory, "] span[ptile=", span.ptile, "]' -u aimin.yan@med.miami.edu",
                             sep = " ")
              },
              general = {
                cmd0 = paste(Wall.time, "-n", cores, "-q general -R 'rusage[mem=",
                             Memory, "] span[ptile=", span.ptile, "]' -u aimin.yan@med.miami.edu",
                             sep = " ")              
              }
      )
      
      job.name = paste0("bamSort.", u)
      cmd1 = paste0("bsub -P bbc -J \"", job.name, paste0("\" -o %J.",
                                                          job.name, ".log "), paste0("-e %J.", job.name, ".err -W"))
      
      
      cmd2=paste("samtools sort",re.out[[u]],file.path(temp3,paste0(file_name, "_sorted")),sep=" ")
      cmd3 = paste(cmd1, cmd0, cmd2, sep = " ")
    } else
    {
      cmd3 = paste("samtools sort",re.out[[u]],file.path(temp3,paste0(file_name, "_sorted")),sep=" ")
    }
    
    cmd <- cmd3
    
    system(cmd)
    
    cmd
    
  },m.id,job.option,Wall.time,cores,Memory,span.ptile,re.out,temp3)
  
  
  cmd.2 <- lapply(1:length(re.out),function(u,m.id,job.option,Wall.time, cores, Memory,
                                            span.ptile,re.out,temp3){
    
    file_name = file_path_sans_ext(basename(re.out[[u]]))
    
    # file_name <- paste0(path_name2,'-',file_name)
    #u <- 3
    if (m.id == 1)
    {
      
      #job.option <- match.arg(job.option)
      
      switch (job.option,
              parallel = {
                cmd0 = paste(Wall.time, "-n", cores, "-q parallel -R 'rusage[mem=",
                             Memory, "] span[ptile=", span.ptile, "]' -u aimin.yan@med.miami.edu",
                             sep = " ")        
              },
              bigmem = {
                cmd0 = paste(Wall.time, "-n", cores, "-q bigmem -R 'rusage[mem=",
                             Memory, "] span[ptile=", span.ptile, "]' -u aimin.yan@med.miami.edu",
                             sep = " ")
              },
              general = {
                cmd0 = paste(Wall.time, "-n", cores, "-q general -R 'rusage[mem=",
                             Memory, "] span[ptile=", span.ptile, "]' -u aimin.yan@med.miami.edu",
                             sep = " ")              
              }
      )
      
      job.name = paste0("bamIndex.", u)
      wait.job.name = paste0("bamSort.", u)
      
      if(wait == TRUE){
        cmd1 = paste0("bsub -w \"done(\"", wait.job.name, "\")\"", " -P bbc -J \"",
                      job.name, paste0("\" -o %J.", job.name, ".log "), paste0("-e %J.",
                                                                               job.name, ".err -W"))
      }else{
        cmd1 = paste0("bsub -P bbc -J \"",job.name, paste0("\" -o %J.", job.name, ".log "), paste0("-e %J.",
                                                                                                   job.name, ".err -W"))
      }
      cmd2=paste("samtools index",file.path(temp3,paste0(file_name, "_sorted.bam")),sep=" ")
      cmd3 = paste(cmd1, cmd0, cmd2, sep = " ")
    } else
    {
      cmd3 = paste("samtools index",file.path(temp3,paste0(file_name, "_sorted.bam")),sep=" ")
    }
    
    cmd <- cmd3
    
    system(cmd)
    
    cmd
    
  },m.id,job.option,Wall.time,cores,Memory,span.ptile,re.out,temp3)
  
  
  cmd.3 <- lapply(1:length(re.out),function(u,m.id,job.option,Wall.time, cores, Memory,
                                            span.ptile,re.out,temp3){
    
    file_name = file_path_sans_ext(basename(re.out[[u]]))
    
    # file_name <- paste0(path_name2,'-',file_name)
    #u <- 3
    if (m.id == 1)
    {
      #job.option <- match.arg(job.option)
      
      switch (job.option,
              parallel = {
                cmd0 = paste(Wall.time, "-n", cores, "-q parallel -R 'rusage[mem=",
                             Memory, "] span[ptile=", span.ptile, "]' -u aimin.yan@med.miami.edu",
                             sep = " ")        
              },
              bigmem = {
                cmd0 = paste(Wall.time, "-n", cores, "-q bigmem -R 'rusage[mem=",
                             Memory, "] span[ptile=", span.ptile, "]' -u aimin.yan@med.miami.edu",
                             sep = " ")
              },
              general = {
                cmd0 = paste(Wall.time, "-n", cores, "-q general -R 'rusage[mem=",
                             Memory, "] span[ptile=", span.ptile, "]' -u aimin.yan@med.miami.edu",
                             sep = " ")              
              }
      )
      
      job.name = paste0("bamPlot.", u)
      wait.job.name = paste0("bamIndex.", u)
      cmd1 = paste0("bsub -w \"done(\"", wait.job.name, "\")\"", " -P bbc -J \"",
                    job.name, paste0("\" -o %J.", job.name, ".log "), paste0("-e %J.",
                                                                             job.name, ".err -W"))
      cmd5="ngs.plot.r -G hg19 -R tss -C"
      cmd6="-O"
      cmd7="-L 4000"
      cmd8=paste(cmd5,file.path(temp3,paste0(file_name,"_sorted.bam")),cmd6,file.path(temp3,paste0(file_name,"_sorted")),cmd7,sep=" ")
      
      
      cmd3 = paste(cmd1, cmd0, cmd8, sep = " ")
    } else
    {
      
      
      "ngs.plot.r -G hg19 -R tss -C config.txt -O 2017-03-02-01_S5_R1.tss -T 2017-03-02-01_S5_R1 -L 4000 -RR 1"
      "ngs.plot.r -G hg19 -R tss -C config.txt -O 2017-03-02-01_S5_R1.tss -T 2017-03-02-01_S5_R1 -L 4000 -RR 1 -GO max -SC global"
      
      cmd5="ngs.plot.r -G hg19 -R tss -C"
      cmd6="-O"
      cmd7="-L 4000"
      cmd8=paste(cmd5,file.path(temp3,paste0(file_name,"_sorted.bam")),cmd6,file.path(temp3,paste0(file_name,"_sorted")),cmd7,sep=" ")
      
      cmd3 = cmd8
    }
    
    cmd <- cmd3
    
    system(cmd)
    
    cmd
    
  },m.id,job.option,Wall.time,cores,Memory,span.ptile,re.out,temp3)
  
  #re <- list(cmd.1 = cmd.1,cmd.2=cmd.2,cmd.3=cmd.3)
  
  # return(re)
}

#' DrawVenn
#'
#' @param devtools
#' @param install_github
#' @param Vennerable
#' @param venneuler
#' @param VennDiagram
#'
#' @return
#' @export
#'
#' @examples
#'
#' out.dir.name="/media/H_driver/2016/Yang/Results/"
#' DrawVenn(out.dir.name)
#'
DrawVenn <- function(out.dir.name) {
  
  devtools::install_github("js229/Vennerable")
  library(Vennerable)
  
  BB=c(1:6153,6154:25736)
  
  CC=c(6154:25736,25737:29579)
  
  AA=c(6026:6153,6154:8462,25737:25776,29650:29719)
  
  Vstem <- Venn(list(ASXL1=AA,SMC1A=BB,RAD21=CC))
  
  Indicator <- (data.matrix(do.call(expand.grid, lapply(seq(1, 
                                                            length = 2), function(x) {
                                                              c(0, 1)
                                                            })))) == 1
  
  
  
  colnames(Indicator) <- SetNames
  
  class(Vstem@IndicatorWeight) <- 
  
  #SetLabels <- VennGetSetLabels(Vstem)
  
  Cstem3 <- compute.Venn(Vstem,doWeights=TRUE)
  
  SetLabels <- VennGetSetLabels(Cstem3)
  FaceLabels <- VennGetFaceLabels(Cstem3)
  
  FaceLabels[FaceLabels$FaceName=="010","x"] <- 10
  FaceLabels[FaceLabels$FaceName=="010","y"] <- 3
  
  SetLabels[SetLabels$Label=="SMC1A","hjust"] <- "right"
  SetLabels[SetLabels$Label=="SMC1A","x"] <- 110
  SetLabels[SetLabels$Label=="SMC1A","y"] <- 58.78653
  
  SetLabels[SetLabels$Label=="RAD21","x"] <- -80
  SetLabels[SetLabels$Label=="RAD21","y"] <- 58.78653
  
  SetLabels[SetLabels$Label=="ASXL1","y"] <-80
  
  Cstem3 <- VennSetSetLabels(Cstem3,SetLabels)
  
  Cstem3<-VennSetFaceLabels(Cstem3,FaceLabels)
  
  Cstem3@FaceLabels[which(Cstem3@FaceLabels$FaceName=="010"),]$x<-90
  Cstem3@FaceLabels[which(Cstem3@FaceLabels$FaceName=="010"),]$y<--35
  
  Cstem3@FaceLabels[which(Cstem3@FaceLabels$FaceName=="001"),]$x<--90
  Cstem3@FaceLabels[which(Cstem3@FaceLabels$FaceName=="001"),]$y<--35
  
  Cstem3@FaceLabels[which(Cstem3@FaceLabels$FaceName=="101"),]$x<--18
  #Cstem3@FaceLabels[which(Cstem3@FaceLabels$FaceName=="101"),]$y<--
  Cstem3@FaceLabels[which(Cstem3@FaceLabels$FaceName=="101"),]$vjust<-"center"
  
  Cstem3@FaceLabels[which(Cstem3@FaceLabels$FaceName=="110"),]$x<-19
  Cstem3@FaceLabels[which(Cstem3@FaceLabels$FaceName=="110"),]$y<-49.86767
  
  c.3.set=sum(Cstem3@IndicatorWeight[,4])
  Cstem3@IndicatorWeight[which(row.names(Cstem3@IndicatorWeight)=="000"),4]=35000-c.3.set
  
  Cstem3@FaceLabels[which(Cstem3@FaceLabels$FaceName=="DarkMatter"),]$y<--130
  
  temp<-ol$venn_cnt
  
  index.A=grep(colnames(Cstem3@IndicatorWeight)[1],colnames(temp))
  index.B=grep("SMC1",colnames(temp))
  index.C=grep("Rad21",colnames(temp))
  
  temp2<-as.data.frame(temp[,c(index.A,index.B,index.C,4)])
  colnames(temp2)[1:3]=colnames(Cstem3@IndicatorWeight)[1:3]
  
  row.names(temp2)<-with(temp2,paste0(temp2$ASXL1,temp2$SMC1A,temp2$RAD21))
  
  row.names(temp2)<-row.names(Cstem3@IndicatorWeight)
  
  temp3<-merge(Cstem3@IndicatorWeight,temp2,by=0,sort=FALSE)
  
  Cstem3@IndicatorWeight[,4]<-temp3$Counts
  
  c.3.set=sum(Cstem3@IndicatorWeight[,4])
  Cstem3@IndicatorWeight[which(row.names(Cstem3@IndicatorWeight)=="000"),4]=35000-c.3.set
  
  library(grid)
  
  grid.newpage()
  
  pdf(paste0(out.dir.name,"venn.pdf"))
  plot(Cstem3)
  dev.off()
  
}

#'  GetResultsFromDiffBind
#'  
#'  Use DiffBind to process peak profiles
#'
#' @param mcf7
#' @param Mergereplicates
#' @param output.file.dir 
#'
#' @return
#' @export
#'
#' @examples
#' output.file.dir="/Volumes/Bioinformatics$/2017/DannyNewData/BindDiff_7_21_2017_2"
#' mcf7=resmcf
#' re<-GetResultsFromDiffBind(mcf7,"yes",output.file.dir)
#' re<-GetResultsFromDiffBind(mcf2,"yes",output.file.dir)
#' 
GetResultsFromDiffBind<-function(mcf7,Mergereplicates=c("yes","no"),output.file.dir){
  
  if (!dir.exists(output.file.dir))
  {
    dir.create(output.file.dir, recursive = TRUE)
  }
  
  temp<-mcf7
  
  sampID.v<-colnames(temp$class) 
  sampID.v.2<-unlist(lapply(1:length(sampID.v),function(u,sampID.v){
    
    x=sampID.v[u]
    
    y=x
    y
  },sampID.v))
  
  
  colnames(temp$class)<-sampID.v.2
  temp$class[1,]<-sampID.v.2
  
  #temp<-dba(temp,mask =c(1,2,13,14,11,12,15,16))
  
  #Merge the replicates of each set 
  if(Mergereplicates=="yes"){
    temp2<-dba.peakset(temp,consensus = -DBA_REPLICATE)
    #temp22<-dba(temp2,mask =c(9,16,17:23))
    print(temp2)
    ##need to set the flexiable number identified from data sets
    t <- temp2
    temp22<-dba(t,mask =which(t$class[which(rownames(t$class)=="Replicate"),]=="1-2"))
  }else{
    temp22<-temp
  }
  
  print(temp22)
  
  Tissue1<-temp22$class[2,]
  Tissue2<-unique(temp22$class[2,])
  
  TF<-unique(temp22$class[3,])
  TF.n<-length(TF)
  
  temp3=file.path(output.file.dir,"Venn")
  
  if(!dir.exists(temp3))
  {
    dir.create(temp3)
  }
  
  for(i in 1:length(Tissue2)){
    po<-which(Tissue1 %in% Tissue2[i])
    
    print(po)
    
    if(length(po)==2)
    {
      png(file.path(temp3,paste0(paste0(colnames(temp22$class)[po],collapse = "-vs-"),".png")))
      dba.plotVenn(temp22,mask=po,main=Tissue2[i])
      dev.off()
      
    }else if(length(po)==4){
      po1<-po[c(1,3)]
      po2<-po[c(2,4)]
      
      png(file.path(temp3,paste0(paste0(colnames(temp22$class)[po1],collapse = "-vs-"),".png")))
      dba.plotVenn(temp22,mask=po1,main=Tissue2[i])
      dev.off()
      
      png(file.path(temp3,paste0(paste0(colnames(temp22$class)[po2],collapse = "-vs-"),".png")))
      dba.plotVenn(temp22,mask=po2,main=Tissue2[i])
      dev.off()
      
    }else
    {
      cat(paste0("For ",Tissue2[i],": Only one peak profile for ",TF.n," TFs\n"))
    }
  }
  
  p.common<-lapply(1:length(Tissue2),function(u,Tissue1,Tissue2,temp22){
    
    po<-which(Tissue1 %in% Tissue2[u])
    
    if(length(po)==2)
    {
      common.peaks<-dba.overlap(temp22,mask=po)
      y<-common.peaks$inAll
    }else if(length(po)==4){
      po1<-po[c(1,3)]
      po2<-po[c(2,4)]
      
      common.peaks.1<-dba.overlap(temp22,mask=po1)
      y1<-common.peaks.1$inAll
      
      common.peaks.2<-dba.overlap(temp22,mask=po2)
      y2<-common.peaks.2$inAll
      
      y<-list(y1=y1,y2=y2)
      
      names(y)[1]<-paste0(colnames(temp22$class)[po1],collapse = "-vs-")
      names(y)[2]<-paste0(colnames(temp22$class)[po2],collapse = "-vs-")
      
    }else{
      y<-NULL}
    y
  },Tissue1,Tissue2,temp22)
  
  names(p.common)<-Tissue2
  
  p.common<-unlist(p.common,recursive = F)
  
  p.common<-p.common[lapply(p.common,length) > 0]
  
  if(!dir.exists(file.path(output.file.dir,"common_peaks_bed")))
  {
    dir.create(file.path(output.file.dir,"common_peaks_bed"))
    dir.create(file.path(output.file.dir,"common_peaks_bed","ucsc"))
    dir.create(file.path(output.file.dir,"common_peaks_bed","igv"))
  }
  
  #output common peaks to bed files
  
  lapply(1:length(p.common),function(u,p.common,output.file.dir){
    x=p.common[[u]]
    
    x_name=names(p.common)[u]
    
    df <- data.frame(seqnames=seqnames(x),
                     #starts=start(x)-1,
                     starts=start(x),
                     ends=end(x),
                     names=c(rep(".", length(x))),
                     scores=elementMetadata(x)[,1],
                     strands=strand(x))
    
    #assign strand
    df.str <- data.frame(seqnames=seqnames(x),
                         #starts=start(x)-1,
                         starts=start(x),
                         ends=end(x),
                         names=c(rep(".", length(x))),
                         scores=elementMetadata(x)[,1],
                         strands=c(rep(".", length(x))))
    
    df.str.1<-df.str[-grep("random",df.str$seqnames),]
    
    df.str.2<-df.str.1
    
    df.str.3<-df.str.2[-grep("chrUn",df.str.2$seqnames),]
    
    write.table(df,file=file.path(output.file.dir,"common_peaks_bed",paste0(x_name,"_cp_with_header.bed")),
                col.names=TRUE,row.names = FALSE,quote=FALSE,sep="\t")
    
    write.table(df,file=file.path(output.file.dir,"common_peaks_bed","ucsc",paste0(x_name,"_4_ucsc.bed")),
                col.names=FALSE,row.names = FALSE,quote=FALSE,sep="\t")
    
    write.table(df,file=file.path(output.file.dir,"common_peaks_bed",paste0(x_name,"_common_peaks.bed")),
                col.names=FALSE,row.names = FALSE,quote=FALSE,sep="\t")
    
    write.table(df,file=file.path(output.file.dir,"common_peaks_bed","igv",paste0(x_name,"_4_igv.bed")),
                col.names=FALSE,row.names = FALSE,quote=FALSE,sep="\t")
    
  },p.common,output.file.dir)
  
  AnntationUsingChipSeeker(file.path(output.file.dir,"common_peaks_bed","igv"),"bed",file.path(output.file.dir,"common_peaks_bed")
                           ,txdb="hg19",DD=5000,distanceToTSS_cutoff=10000)
  
  return(p.common)
}

#' GetSampleInfo
#'
#' GetSampleInfo 
#' 
#' @param input.sample.file
#' @param input.bam.file 
#'  
#' 
#' @return
#' @export
#'
#' @examples
#' 
#' input.sample.file <- '/Volumes/Bioinformatics$/2017/DannyNewData/SampleID_INFO_ChIP_new_Danny.csv'
#' 
#' input.bam.file <- '/Volumes/Bioinformatics$/2017/DannyNewData/sorted_bam_files.txt'
#' 
#' input.bam.file <- '/Volumes/Bioinformatics$/2017/DannyNewData/NewRe2Danny/sorted_bam_files_2.txt'
#' 
#' re <- GetSampleInfo(input.sample.file,input.bam.file)
#' 
#' 
GetSampleInfo <- function(input.sample.file, input.bam.file)
{
  
  re <- read.csv(input.sample.file)
  re.c <- colnames(re)
  
  bam.file <- read.table(input.bam.file, header = FALSE)
  
  yy <- apply(bam.file, 1, function(u)
  {
    y <- basename(u)
    pos <- regexpr("\\.", y)
    pos <- pos - 1
    y <- substr(y, 1, pos)
    y
  })
  
  bam.file.sample.name <- cbind(bam.file, yy)
  
  colnames(bam.file.sample.name) <- c("file.name", "ID")
  
  Ab.type <- unique(as.character(re[, colnames(re) == re.c[3]]))
  Cell.type <- unique(as.character(re[, colnames(re) == re.c[2]]))
  
  re1 <- cbind(re[, 1:3], paste0(re[, 2], "_", re[, 3]))
  
  re11 <- merge(re1, bam.file.sample.name, by = "ID")
  
  colnames(re11)[4] = "Cell_TF"
  
  chipseqCell <- function(n, es)
  {
    value <- list(name = n, es = es)
    attr(value, "class") <- "ChipSeqCell"
    value
  }
  
  cellType <- unique(as.character(re1[re1$Type_TF == "Input", ]$Type_Cell))
  
  
  yy <- lapply(cellType, function(u, re11)
  {
    n <- nchar(u)
    index1 <- which(substr(re11$Type_Cell, 1, n) == u)
    z <- re11[index1, ]
    index2 <- which(nchar(as.character(z$Type_Cell)) <= n + 2)
    
    z2 <- z[index2, ]
    
    cc <- chipseqCell(u, z2)
    
    cc
  }, re11)
  
  colnames(re1)[4] = "Cell_TF"
  
  Cell.Ab.type <- unique(as.character(re1[, 4]))
  
  re2 <- lapply(1:length(Cell.Ab.type), function(u, Cell.Ab.type, re1)
  {
    
    x = Cell.Ab.type[u]
    z = re1
    
    ZZ <- as.character(z[which(z[, 4] == x), 1])
    ZZ
    
  }, Cell.Ab.type, re1)
  
  names(re2) <- Cell.Ab.type
  
  re21 <- re2[lapply(re2, length) == 1]
  
  re22 <- re2[lapply(re2, length) == 2]
  
  re3 <- list(re1 = re1, re2 = re2, re21 = re21, re22 = re22, re11 = re11, 
              y = yy)
  
  return(re3)
  
}

#' peakcallwithinput
#' 
#' @export
#' @example 
#' 
#' genome='Hs'
#' re <- peakcallwithinput(input.sample.file,input.bam.file,genome,output.dir,peakcaller)
#' 
#' 
peakcallwithinput <- function(input.sample.file, input.bam.file, genome = c("Hs", 
                                                                            "hs", "HS", "hS"), output.dir, peakcaller = c("macs14", "macs2"), peakPvalue)
{
  
  re <- GetSampleInfo(input.sample.file, input.bam.file)
  
  cellInfo <- re$y
  
  output.dir.name = dirname(input.sample.file)
  
  temp3 = file.path(output.dir.name, output.dir)
  
  if (!dir.exists(temp3))
  {
    dir.create(temp3)
  }
  
  peakcaller <- match.arg(peakcaller)
  
  genome <- match.arg(genome)
  cmd10 <- paste("-f BAM", "-g", genome, "-n", sep = " ")
  
  switch(peakcaller, macs2 = {
    
    PATH1 = Sys.getenv("PATH")
    
    macs2_Lib = file.path("/nethome/axy148/NGS_tools/MACS/bin/")
    
    Sys.setenv(PATH = paste0(macs2_Lib, ":", PATH1))
    
    cmd1 <- Sys.which("macs2")[[1]]
    
    cat(cmd1, "\n")
    
    cmd9 = paste(cmd1, "callpeak -t", sep = " ")
    
    cmd11 <- paste("-p", peakPvalue, sep = " ")
    
  }, {
    cmd9 = "macs14 -t "
    
    cmd11 <- paste("-p", peakPvalue, sep = " ")
    
  })
  
  cellInfo.run <- lapply(1:length(cellInfo), function(u, cellInfo, temp3)
  {
    
    x.name = cellInfo[[u]]$name
    
    es <- cellInfo[[u]]$es
    
    x.input <- es[es$Type_TF == "Input", ]$file.name
    
    x.sample <- es[es$Type_TF != "Input", ]
    
    x.run <- apply(x.sample, 1, function(x)
    {
      
      y <- x
      
      ID <- y[1]
      Type_Cell <- y[2]
      Type_TF <- y[3]
      Cell_TF <- y[4]
      file.name <- y[5]
      xx <- file.name
      xx.name = paste(ID, gsub(" ", "-", Type_Cell), Type_TF, sep = "-")
      
      cmd12 = paste(cmd9, xx, "-c", x.input, cmd10, file.path(temp3, paste0(xx.name, 
                                                                            "_hs_1.00e-05_", peakcaller)), cmd11, sep = " ")
      
      cmd12
      
    })
    
    x.run
    
  }, cellInfo, temp3)
  
  
  names(cellInfo.run) = unlist(lapply(cellInfo, function(u)
  {
    u$name
  }))
  
  zzz <- unlist(cellInfo.run)
  
  lapply(1:length(zzz), function(u, zzz)
  {
    
    cat(as.character(zzz[u][[1]]), "\n")
    cat("\n")
    
    system(as.character(zzz[u][[1]]))
    
  }, zzz)
  
  
  re <- list(cellInforun = cellInfo.run, zzz = zzz)
  
  AnntationUsingChipSeeker(temp3, "peaks.bed", temp3, DD = 5000)
  
  return(re)
  
}

select.sample <- c("MDA MB 231-DD-1_cJun", "MDA MB 231-1_cJun", "1833-1_cJun")
output.config.dir <- "~/"

configAndMultiplot <- function(res, select.sample, output.config.dir)
{
  
  if (!dir.exists(output.config.dir))
  {
    dir.create(output.config.dir)
  }
  
  x <- res$re11
  
  xx <- x[which(x$Cell_TF %in% select.sample), ]
  
  xxx <- cbind.data.frame(xx$file.name, rep("-1", dim(xx)[1]), gsub(" ", "-", 
                                                                    xx$Cell_TF))
  
  config.sample.name <- paste(gsub(" ", "-", select.sample), collapse = "-and-")
  
  config.file.name <- file.path(output.config.dir, paste0(config.sample.name, 
                                                          "-config.txt"))
  
  write.table(xxx, file = config.file.name, col.names = FALSE, row.names = FALSE, 
              quote = FALSE, sep = "\t")
  
  cmd0 <- "ngs.plot.r -G hg19 -R tss -C"
  
  cmd1 <- "-O"
  
  cmd2 <- "-L 4000 -RR 1 -CD 1 -CO \"blue\""
  
  output.results <- file.path(output.config.dir, paste0(config.sample.name, 
                                                        "results"))
  
  cmd <- paste(cmd0, config.file.name, cmd1, output.results, cmd2, sep = " ")
  
  system(cmd)
  
  return(config.file.name)
  
}

#' input.sample.file <- '/Volumes/Bioinformatics$/2017/DannyNewData/SampleID_INFO_ChIP_new_Danny.csv'
#' 
#' input.bam.file <- '/Volumes/Bioinformatics$/2017/DannyNewData/NewRe2Danny/sorted_bam_files_2.txt'
#' 
#' re <- ChipSeq:::matchBamInputGene(input.sample.file,input.bam.file,'$HOME/all_common_gene_unique.txt','$HOME/NgsConfigFile')
#' 
#' R -e 'libraray(ChipSeq);re <- ChipSeq:::matchBamInputGene('/scratch/projects/bbc/aiminy_project/DannyNewData2/SampleID_INFO_ChIP_new_Danny.csv','/scratch/projects/bbc/aiminy_project/DannyNewData2/sorted_bam_files_2.txt','$HOME/all_common_gene_unique.txt','$HOME/NgsConfigFile',,job.option = 'parallel')'
#' 
#' bsub -P bbc -J 'bamPlot' -o %J.bamPlot.log -e %J.bamPlot.err -W 72:00 -n 32 -q parallel -R 'rusage[mem= 16000 ] span[ptile= 16 ]' -u aimin.yan@med.miami.edu R -e 'library(ChipSeq);re <- ChipSeq:::matchBamInputGene('/scratch/projects/bbc/aiminy_project/DannyNewData2/SampleID_INFO_ChIP_new_Danny.csv','/scratch/projects/bbc/aiminy_project/DannyNewData2/sorted_bam_files_2.txt','~/all_common_gene_unique.txt','NgsConfigFile')'
#' 
#' bsub -P bbc -J "bamPlot" -o %J.bamPlot.log -e %J.bamPlot.err -W 72:00 -n 32 -q parallel -R 'rusage[mem= 16000 ] span[ptile= 16 ]' -u aimin.yan@med.miami.edu R -e 'library(ChipSeq);re <- ChipSeq:::matchBamInputGene("/scratch/projects/bbc/aiminy_project/DannyNewData2/SampleID_INFO_ChIP_new_Danny.csv","/scratch/projects/bbc/aiminy_project/DannyNewData2/sorted_bam_files_2.txt","~/all_common_gene_unique.txt","NgsConfigFile3",add.input= "yes")'
matchBamInputGene <- function(input.sample.file, input.bam.file, input.gene.list, 
                              output.dir, ngs.para = c("hg19", 4000, 1, 1, "total"), add.input = NULL)
{
  
  re <- GetSampleInfo(input.sample.file, input.bam.file)
  
  cellInfo <- re$y
  
  # output.dir.name = dirname(input.sample.file)
  
  if (!dir.exists(output.dir))
  {
    dir.create(output.dir, recursive = TRUE)
  }
  
  temp3 = output.dir
  
  # cmd9 = 'ngs.plot.r -G' cmd10 = '-R' cmd11 = '-C' cmd12 = '-O' cmd13 = '-T'
  # cmd14 = '-L' cmd15 = '-RR' cmd16 = '-CD' cmd17= '-GO'
  
  
  
  cellInfo.run <- lapply(1:length(cellInfo), function(u, add.input, cellInfo, 
                                                      temp3)
  {
    
    x.name = cellInfo[[u]]$name
    
    es <- cellInfo[[u]]$es
    
    x.input <- es[es$Type_TF == "Input", ]$file.name
    
    x.sample <- es[es$Type_TF != "Input", ]
    
    x.run <- apply(x.sample, 1, function(x, add.input, temp3)
    {
      
      y <- x
      
      ID <- y[1]
      Type_Cell <- y[2]
      Type_TF <- y[3]
      Cell_TF <- y[4]
      file.name <- y[5]
      xx <- file.name
      xx.name = paste(ID, gsub(" ", "-", Type_Cell), Type_TF, sep = "-")
      
      
      if (!is.null(add.input))
      {
        cmd12 = paste(paste0(xx, ":", x.input), input.gene.list, paste0(gsub(" ", 
                                                                             "-", Cell_TF), "_cJunAndp27"), sep = "\t")
      } else
      {
        cmd12 = paste(xx, input.gene.list, paste0(gsub(" ", 
                                                       "-", Cell_TF), "_cJunAndp27"), sep = "\t")
      }
      cat(cmd12, file = file.path(temp3, paste0(xx.name, "_config_cJunAndp27.txt")), 
          sep = "\t")
    }, add.input, temp3)
    
    x.run
    
  }, add.input, cellInfo, temp3)
  
  
  # dir.name=temp3 dir.name=reformatPath(dir.name)
  
  file.name = file.path(temp3, dir(temp3, recursive = TRUE))
  
  file.name.2 <- as.list(file.name)
  
  
  # names(file.name.2) = unlist(lapply(file.name.2, function(u) { u$name }))
  
  zzz <- unlist(file.name.2)
  
  lapply(1:length(zzz), function(u, zzz)
  {
    
    dir.name = dirname(zzz[u][[1]])
    file_name = file_path_sans_ext(basename(zzz[u][[1]]))
    
    cmd = paste("ngs.plot.r -G hg19 -R tss -C", zzz[u][[1]], "-O", file.path(dir.name, 
                                                                             paste0(file_name, ".tss")), "-T", file_name, "-L 4000 -RR 1 -CD 1 -GO total", 
                sep = " ")
    
    
    
    # system(as.character(zzz[u][[1]])) job.name = paste0('bamPlot.', u)
    # cmd.pegasus = usePegasus(job.option, Wall.time = '72:00',cores = 32,Memory
    # = 16000,span.ptile = 16,job.name) cmd2 = paste(cmd.pegasus,cmd,sep = ' ')
    
    cmd2 = cmd
    cat(cmd2, "\n")
    cat("\n")
    system(cmd2)
  }, zzz)
  
  
  re <- list(cellInforun = cellInfo.run, zzz = zzz)
  
  # AnntationUsingChipSeeker(temp3, 'peaks.bed', temp3, DD = 5000)
  
  return(re)
  
}

# /deepTools-1.5/bin/bamCompare --bamfile1 ChIP.bam --bamfile2 Input.bam \
# --binSize 25 --fragmentLength 200 --missingDataAsZero no \
# --ratio log2 --scaleFactorsMethod SES -o log2ratio_ChIP_vs_Input.bw

#' bsub -P bbc -J "bamCompare" -o %J.bamCompare.log -e %J.Compare.err -W 72:00 -n 32 -q parallel -R 'rusage[mem= 16000 ] span[ptile= 16 ]' -u aimin.yan@med.miami.edu R -e 'library(ChipSeq);re <- ChipSeq:::useBamCompare("/scratch/projects/bbc/aiminy_project/DannyNewData2/SampleID_INFO_ChIP_new_Danny.csv","/scratch/projects/bbc/aiminy_project/DannyNewData2/sorted_bam_files_2.txt","~/BamCompare")'

useBamCompare <- function(input.sample.file,input.bam.file,output.dir)
{
  
  re <- GetSampleInfo(input.sample.file, input.bam.file)
  
  cellInfo <- re$y
  
  # output.dir.name = dirname(input.sample.file)
  
  if (!dir.exists(output.dir))
  {
    dir.create(output.dir, recursive = TRUE)
  }
  
  temp3 = output.dir
  
  # cmd9 = 'ngs.plot.r -G' cmd10 = '-R' cmd11 = '-C' cmd12 = '-O' cmd13 = '-T'
  # cmd14 = '-L' cmd15 = '-RR' cmd16 = '-CD' cmd17= '-GO'
  
  
  
  cellInfo.run <- lapply(1:length(cellInfo), function(u,cellInfo, 
                                                      temp3)
  {
    
    x.name = cellInfo[[u]]$name
    
    es <- cellInfo[[u]]$es
    
    x.input <- es[es$Type_TF == "Input", ]$file.name
    
    x.sample <- es[es$Type_TF != "Input", ]
    
    #print(x.sample)
    #print(x.input)
    
    x.run <- apply(x.sample, 1, function(x, x.input, temp3)
    {
      
      y <- x
      
      ID <- y[1]
      Type_Cell <- y[2]
      Type_TF <- y[3]
      Cell_TF <- y[4]
      file.name <- y[5]
      xx <- file.name
      xx.name = paste(ID, gsub(" ", "-", Type_Cell), Type_TF, sep = "-")
      
      # ~/python/Python-2.7.11/python  ~/NGS_tools/deepTools/bin/bamCompare  --bamfile1 /scratch/projects/bbc/aiminy_project/DannyNewNgsPlot/2017-03-02-03_S11_R1.marked_sorted.bam --bamfile2 /scratch/projects/bbc/aiminy_project/DannyNewNgsPlot/2017-03-02-17_S1_R1.marked_sorted.bam --binSize 25 --ratio log2 -o ~/BamCompare/log2ratio_2017-03-02-03_S11_R1.marked_sorted.bam_vs_2017-03-02-17_S1_R1.marked_sorted.bam.bw 
      
      cmd1 <- paste("~/python/Python-2.7.11/python  ~/NGS_tools/deepTools/bin/bamCompare --bamfile1",xx,"--bamfile2",x.input,sep=" ")
      cmd2 <- "--binSize 25"
      cmd3 <- "--ratio log2 -o"
      
      cmd4 <- paste(cmd1,cmd2,cmd3,sep=" ") 
      
      cmd5 <- file.path(output.dir,paste0("log2ratio_",basename(as.character(xx)),"_vs_",basename(as.character(x.input)),".bw"))
      
      
      cmd6 <- paste(cmd4,cmd5,sep=" ") 
      
      cmd6
      
      #cat(cmd6, "\n")
      #cat("\n")
    }, x.input, temp3)
    
    x.run
    
  }, cellInfo,temp3)
  
  
  names(cellInfo.run) = unlist(lapply(cellInfo, function(u)
  {
    u$name
  }))
  
  sysyem("module unload python/2.7.3")
  
  zzz <- unlist(cellInfo.run)
  
  lapply(1:length(zzz), function(u, zzz)
  {
    
    cat(as.character(zzz[u][[1]]), "\n")
    cat("\n")
    
    system(as.character(zzz[u][[1]]))
    
  }, zzz)
  # # dir.name=temp3 dir.name=reformatPath(dir.name)
  # 
  # file.name = file.path(temp3, dir(temp3, recursive = TRUE))
  # 
  # file.name.2 <- as.list(file.name)
  # 
  # 
  # # names(file.name.2) = unlist(lapply(file.name.2, function(u) { u$name }))
  # 
  # zzz <- unlist(file.name.2)
  # 
  # lapply(1:length(zzz), function(u, zzz)
  # {
  #   
  #   dir.name = dirname(zzz[u][[1]])
  #   file_name = file_path_sans_ext(basename(zzz[u][[1]]))
  #   
  #   cmd = paste("ngs.plot.r -G hg19 -R tss -C", zzz[u][[1]], "-O", file.path(dir.name, 
  #                                                                            paste0(file_name, ".tss")), "-T", file_name, "-L 4000 -RR 1 -CD 1 -GO total", 
  #               sep = " ")
  #   
  #   
  #   
  #   # system(as.character(zzz[u][[1]])) job.name = paste0('bamPlot.', u)
  #   # cmd.pegasus = usePegasus(job.option, Wall.time = '72:00',cores = 32,Memory
  #   # = 16000,span.ptile = 16,job.name) cmd2 = paste(cmd.pegasus,cmd,sep = ' ')
  #   
  #   cmd2 = cmd
  #   cat(cmd2, "\n")
  #   cat("\n")
  #   system(cmd2)
  # }, zzz)
  # 
  # 
  # re <- list(cellInforun = cellInfo.run, zzz = zzz)
  # 
  # # AnntationUsingChipSeeker(temp3, 'peaks.bed', temp3, DD = 5000)
  
  #return(re)
  
}

#' x <- ChipSeq:::usePegasus('parallel', Wall.time = '72:00',cores = 32,Memory = 16000,span.ptile = 16,'bamPlot')

usePegasus <- function(job.option = c("general", "parallel", "bigmem"), Wall.time, 
                       cores, Memory, span.ptile, job.name, wait.job.name = NULL)
{
  
  job.option <- match.arg(job.option)
  
  switch(job.option, parallel = {
    cmd0 = paste(Wall.time, "-n", cores, "-q parallel -R 'rusage[mem=", 
                 Memory, "] span[ptile=", span.ptile, "]' -u aimin.yan@med.miami.edu", 
                 sep = " ")
  }, bigmem = {
    cmd0 = paste(Wall.time, "-n", cores, "-q bigmem -R 'rusage[mem=", Memory, 
                 "] span[ptile=", span.ptile, "]' -u aimin.yan@med.miami.edu", sep = " ")
  }, general = {
    cmd0 = paste(Wall.time, "-n", cores, "-q general -R 'rusage[mem=", Memory, 
                 "] span[ptile=", span.ptile, "]' -u aimin.yan@med.miami.edu", sep = " ")
  })
  
  if (!is.null(wait.job.name))
  {
    cmd1 = paste0("bsub -w \"done(\"", wait.job.name, "\")\"", " -P bbc -J \"", 
                  job.name, paste0("\" -o %J.", job.name, ".log "), paste0("-e %J.", 
                                                                           job.name, ".err -W"))
  } else
  {
    cmd1 = paste0("bsub -P bbc -J \"", job.name, paste0("\" -o %J.", job.name, 
                                                        ".log "), paste0("-e %J.", job.name, ".err -W"))
  }
  
  cmd = paste(cmd1, cmd0, sep = " ")
  
  return(cmd)
}

#' IndexBamFile
#'
#' @param input.file.dir 
#' @param input.file.pattern 
#' @param index.file 
#' @param output.file.dir 
#' @param genome 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' input.file.dir="/projects/scratch/bbc/Project/Danny_chip/Alignment/BWA/"
#' input.file.pattern="*.bam"
#' index.file=9
#' output.file.dir="/scratch/projects/bbc/aiminy_project/
#' genome="Hs"
#' 
#' IndexBamFile(input.file.dir,input.file.pattern,index.file,output.file.dir,genome)
#' 
IndexBamFile <- function(input.file.dir,input.file.pattern,index.file,output.file.dir,genome) {
  
  #library(ChIPpeakAnno)
  
  dir.name=input.file.dir
  input.file.pattern=input.file.pattern
  
  dir.name=reformatPath(dir.name)
  output.dir.name=reformatPath(output.file.dir)
  
  #print(output.dir.name)
  temp=Sys.time()
  temp1=gsub(":","-",Sys.time())
  temp2=gsub(" ","-",temp1)
  temp3=paste0(output.dir.name,"AnalysisResults_at_",temp2)
  
  dir.create(temp3)
  
  file.name=paste0(dir.name,dir(dir.name,recursive = TRUE,pattern=input.file.pattern))
  file.name.2<-as.list(file.name)
  
  names(file.name.2)=sapply(strsplit(file.name,split="\\/"),"[[",index.file)
  
  #print(file.name.2)
  
  re.out<-file.name.2
  
  cmd1="samtools index "
  
  lapply(1:length(re.out),function(u,re.out,temp3){
    
    x=re.out[[u]]
    x_name=names(re.out)[u]
    cmd2=paste0(cmd1," ",x," ",paste0(temp3,"/",x_name,"_sorted"))
    
    #print(cmd2)
    system(cmd2)
    
  },re.out,temp3)
  
  #sample.name<-sapply(strsplit(names(file.name.3),split="_peaks_"),"[[",1)
  
  #names(file.name.3)=sample.name
  
  #file.name.4 <-file.name.3[-1]
  
  #re.out<-lapply(file.name.4,function(u){
  #  re=toGRanges(u,format="BED")
  #colnames(re)=c("Count","GeneName")
  #  re
  #})
  
  #head(re.out[[1]])
  
  #   re.out.L<-lapply(re.out,function(u){
  #     re=length(u)
  #     #colnames(re)=c("Count","GeneName")
  #     re
  #   })
  #   
  #   
  # #  cmd="samtools sort"
  # #  input.file="/scratch/projects/bbc/Project/Danny_chip2/Alignment/BWA/""2016-10-14-Hyunho-Yoon-16-1833-PF1502-1-p27_S4_R1.bam" 
  # #  output.file="/scratch/projects/bbc/aiminy_project/Bam_marked_sorted/""2016-10-14-Hyunho-Yoon-16-1833-PF1502-1-p27_S4_R1.bam_sorted"
  #   
  #   if(genome=="Mm"){
  #     
  #     annoData <- toGRanges(EnsDb.Mmusculus.v75, feature="gene")
  #     
  #     ol <- findOverlapsOfPeaks(re.out[c(2,4,1)])
  #     
  #     overlaps<-ol$peaklist$`11_2470IUPUI_WT_BM_SMC1_peaks.bed///13_2470IUPUI_WT_BM_Rad21_peaks.bed///10_WT_BM_ASXL1_peaks.bed`
  #     
  #     binOverFeature(overlaps, annotationData=annoData,
  #                    radius=5000, nbins=20, FUN=length, errFun=0,
  #                    ylab="count",
  #                    main="Distribution of aggregated peak numbers around TSS")
  #     
  #     overlaps.trimmed<-trim(overlaps, use.names=TRUE)
  #     
  #     library(EnsDb.Mmusculus.v79)
  #     dd.GRCm39.mm10<-toGRanges(EnsDb.Mmusculus.v75)
  #     #seqinfo(dd.GRCm39.mm10)
  #     #seqlevels(dd.GRCm39.mm10)
  #     
  #     seqlevels(dd.GRCm39.mm10,force=TRUE) <- c("chr1","chr10","chr11","chr12","chr13",
  #                                               "chr14","chr15","chr16","chr17","chr18","chr19","chr2",
  #                                               "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY")
  #     
  #     seqinfo(overlaps)<-seqinfo(dd.GRCm39.mm10)
  #     #GRCm38/mm10
  #     #dd<-toGRanges(EnsDb.Mmusculus.v79)
  #     #seqinfo(dd)
  #     #library(ensembldb)
  #     #library(GenomeInfoDb)
  #     seqlevelsStyle(overlaps.trimmed) <- seqlevelsStyle(dd.GRCm39.mm10)
  #     overlaps.anno<-annoPeaks(overlaps.trimmed,dd.GRCm39.mm10)
  #     write.table(overlaps.anno,file=paste0(temp3,"/","annotation.txt"),row.names = FALSE,quote=FALSE,sep="\t")
  #   }else if(genome=="Hs"){
  #     
  #     library(EnsDb.Hsapiens.v75)
  #     #annoData<-toGRanges(EnsDb.Hsapiens.v75, feature="gene")
  #     
  #     dd.hs<-toGRanges(EnsDb.Hsapiens.v75)
  #     
  #     print(seqinfo(dd.hs))
  #     print(seqlevels(dd.hs))
  #     
  #     #print(seqlevels(dd.hs)[,1])
  #     
  #     #print(seqlevels(re.out[[1]])[,1])
  #     
  #     # seqlevels(dd.hs,force=TRUE) <- c("chr1","chr10","chr11","chr12","chr13",
  #     #                                           "chr14","chr15","chr16","chr17","chr18","chr19","chr2",
  #     #                                           "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY")
  #     
  #     #temp4=
  #     
  #     re.out.L<-lapply(1:length(re.out),function(u,re.out,dd.hs){
  #       
  #       x=re.out[[u]]
  #       x_name=names(re.out)[u]
  #       
  #       seqlevels(dd.hs,force=TRUE)<-seqinfo(x)@seqnames
  #       #print(seqinfo(re.out.trimmed))
  #       #print(seqlevels(re.out.trimmed))
  #       seqinfo(x)<-seqinfo(dd.hs)
  #       #GRCm38/mm10
  #       #dd<-toGRanges(EnsDb.Mmusculus.v79)
  #       #seqinfo(dd)
  #       #library(ensembldb)
  #       #library(GenomeInfoDb)
  #       seqlevelsStyle(x) <- seqlevelsStyle(dd.hs)
  #       re.out.trimmed<-trim(x, use.names=TRUE)
  #       overlaps.anno<-annoPeaks(re.out.trimmed,dd.hs)
  #       
  #       write.table(overlaps.anno,file=paste0(temp3,"/",x_name,"_annotation.txt"),row.names = FALSE,quote=FALSE,sep="\t")
  #     },re.out,dd.hs)
  #     
  #   }
  
}

#' InstallRequiredPackage
#'
#' @return
#' @export
#'
#' @examples
#' InstallRequiredPackage()
#'
InstallRequiredPackage <- function() {
  
  #install from bioc
  
  source("http://bioconductor.org/biocLite.R")
  biocLite("GenomicRanges")
  biocLite("ChIPseeker")
  biocLite("regioneR")
  biocLite("DiffBind")
  biocLite("BSgenome.Hsapiens.UCSC.hg19")
  biocLite("TxDb.Mmusculus.UCSC.mm9.knownGene")
  biocLite("EnsDb.Hsapiens.v75")
  devtools::install_github("sheffien/simpleCache")
  biocLite("LOLA")
  biocLite("ChIPpeakAnno")
  biocLite("BSgenome.Mmusculus.UCSC.mm9")
  biocLite("motifStack")
  biocLite("BSgenome.Ecoli.NCBI.20080805")
  biocLite("EnsDb.Mmusculus.v75")
  biocLite("EnsDb.Mmusculus.v79")
  biocLite("BSgenome.Mmusculus.UCSC.mm10")
  biocLite("org.Mm.eg.db")
  
  #install from git
  devtools::install_github("nsheff/LOLA")
  
  
  #loading
  library(BiocInstaller)
  library(BiocInstaller)
  library(ChIPseeker)
  library(DiffBind)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(TxDb.Mmusculus.UCSC.mm9.knownGene)
  library(LOLA)
  library("simpleCache")
  library(ChIPpeakAnno)
  library(rtracklayer)
  library(org.Mm.eg.db)
  library(motifStack)
  
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(BSgenome.Mmusculus.UCSC.mm9)
  library(EnsDb.Mmusculus.v75)
  library(EnsDb.Mmusculus.v79)
  
  library(ggplot2)
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
#' LoadRequiredPackage()
#'
LoadRequiredPackage <- function() {
  library(LOLA)
  
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  library(clusterProfiler)
  library(regioneR)
  library(BSgenome.Hsapiens.UCSC.hg19)
}

#' ParserBamFile4NgsPlot
#'
#' @param dir.name
#' @param input.file.pattern
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#' sh ~/Code/BashRunMACS1-4-2_4_Danny_chip_seq3.sh /scratch/projects/bbc/aiminy_project/Bam_sorted/ 7 ".sorted.bam$" /scratch/projects/bbc/aiminy_project/Bam_marked_sorted/
#'
#' dir.name="/scratch/projects/bbc/aiminy_project/Bam_sorted/"
#' index=7
#' input.file.pattern=".sorted.bam$"
#' out.dir.name="/scratch/projects/bbc/aiminy_project/Bam_marked_sorted/"
#' outfile="test.txt"
#'
#' re.from.bed.peng.4.venn<-ParserBamFile4NgsPlot(dir.name,index,input.file.pattern,out.dir.name,outfile)
#'
#' save.image(file=paste0(out.dir.name,"re_save_peng.RData"))
#' savehistory(file=paste0(out.dir.name,"re_save_peng.Rhistory"))
#'
ParserBamFile4NgsPlot<-function(dir.name,index,input.file.pattern,out.dir.name,outfile){
  
  file.name=paste0(dir.name,dir(dir.name,recursive = TRUE,pattern=input.file.pattern))
  file.name.2<-as.list(file.name)
  
  
  #file.name="/scratch/projects/bbc/Project/Danny_chip/Alignment/BWA/2016-10-14-Hyunho-Yoon-02-MDA-MB-231-2-cJun_S5_R1.marked.bam"
  
  #strsplit(file.name,split="\\/")
  
  names(file.name.2)=sapply(strsplit(file.name,split="\\/"),"[[",index)
  
  print(file.name.2)
  
  ab<-names(file.name.2)
  
  #ab<-"2016-10-14-Hyunho-Yoon-13-MDA-MB-231-DD-2-p27_S15_R1.marked.bam"
  
  ab2<-sapply(strsplit(ab,split="\\."),"[[",1)
  
  print(ab2)
  #ab3<-sapply(strsplit(ab2,split="\\-"),"[[",length(strsplit(ab2,split="\\-")[[1]]))
  
  
  re<-cbind(paste0(dir.name,names(file.name.2)),rep(-1,length(names(file.name.2))),ab2)
  
  write.table(re, file =paste0(out.dir.name,outfile), append = FALSE, quote = F, sep = " ",eol = "\n", na = "NA", dec = ".", row.names = F, col.names = F)
  
  
  cmd1="ngs.plot.r -G hg19 -R tss -C"
  cmd2="-O"
  cmd3="-L 4000"
  #cmd3="-L 4000 -RR 1 -CD 1 -CO \\\"blue\\\""
  
  #ngs.plot.r -G hg19 -R tss -C $1 -O $2 -L 4000 -RR 1 -CD 1 -CO "blue"
  
  
  #file.name.3<-file.name.2[-6]
  
  re.out<-lapply(file.name.2,function(u){
    
    cmd4=paste(cmd1,u,cmd2,u,cmd3,sep=" ")
    
    print(cmd4)
    system(cmd4, intern = TRUE, ignore.stderr = TRUE)
    
    #re=read.table(u,header=FALSE)
    #  re<-as.character(re[,1])
    #  #colnames(re)=c("Count","GeneName")
    #  re
  })
  
  #sample.name<-sapply(strsplit(names(re.out),split="_peaks_"),"[[",1)
  
  #names(re.out)=sample.name
  
  # A=re.out[[1]] #10_2470IUPUI_WT_BM_ASXL1
  # B=re.out[[2]] #10_WT_BM_ASXL1
  # C=re.out[[3]] #11_2470IUPUI_WT_BM_SMC1
  # D=re.out[[5]] #13_2470IUPUI_WT_BM_Rad21
  #
  # A.name=sample.name[1] #10_2470IUPUI_WT_BM_ASXL1
  # B.name=sample.name[2] #10_WT_BM_ASXL1
  # C.name=sample.name[3] #11_2470IUPUI_WT_BM_SMC1
  # D.name=sample.name[5] #13_2470IUPUI_WT_BM_Rad21
  #
  # #A <- createRandomRegions(nregions=50, length.mean=5000000, length.sd=3000000)
  # #B <- c(A[1:25], createRandomRegions(nregions=25, length.mean=500000, length.sd=30000))
  # #class(A)
  # #class(B)
  #
  #
  # # TestOverlap <- function(A, B, A.name,B.name,out.dir.name,n_permutations) {
  # #   NumOverLap<-numOverlaps(A, B, count.once=TRUE)
  # #   pt <- overlapPermTest(A=A, B=B, ntimes=n_permutations)
  # #   png(filename=paste0(out.dir.name,A.name,B.name,n_permutations,"_overlap_test.png"))
  # #   plot(pt)
  # #   dev.off()
  # #
  # #   lz <- localZScore(pt=pt, A=A, B=B)
  # #   png(filename=paste0(out.dir.name,A.name,B.name,n_permutations,"_localZScore.png"))
  # #   plot(lz)
  # #   dev.off()
  # # }
  # #
  # # TestOverlap(A,B,A.name,B.name,out.dir.name,100)
  # # TestOverlap(A,C,A.name,C.name,out.dir.name,100)
  # # TestOverlap(A,D,A.name,D.name,out.dir.name,100)
  # #
  # # TestOverlap(B,C,B.name,C.name,out.dir.name,100)
  # # TestOverlap(B,D,B.name,D.name,out.dir.name,100)
  # # TestOverlap(C,D,C.name,D.name,out.dir.name,100)
  # #
  #
  # # temp.name<-strsplit(names(file.name.2),split="\\.")
  #
  # #a = readBed(system.file("extdata", "examples/combined_regions.bed",
  # #                        package="LOLA"))
  #
  #
  # venn.plot <- venn.diagram(
  #   x = re.out[c(1,2)],
  #   filename = paste0(out.dir.name,names(re.out)[1],names(re.out)[2],"_overlap_venn.tiff"),
  #   height = 3000,
  #   width = 3500,
  #   resolution = 1000,
  #   col = "black",
  #   lty = "dotted",
  #   lwd = 1,
  #   fill = c("red","blue"),
  #   alpha = 0.50,
  #   label.col = c(rep("white",3)),
  #   cex = 0.5,
  #   fontfamily = "serif",
  #   fontface = "bold",
  #   cat.col = c("red","blue"),
  #   cat.cex = 0.5,
  #   cat.pos = 0.5,
  #   cat.dist = 0.05,
  #   cat.fontfamily = "serif"
  # )
  #
  # venn.plot <- venn.diagram(
  #   x = re.out[c(1,3,5)],
  #   filename =paste0(out.dir.name,names(re.out)[1],names(re.out)[3],names(re.out)[5],"_overlap_venn.tiff"),
  #   col = "black",
  #   lty = "dotted",
  #   lwd = 2,
  #   fill = c("red", "orange", "blue"),
  #   alpha = 0.50,
  #   label.col = c(rep("white",7)),
  #   cex = 1,
  #   fontfamily = "serif",
  #   fontface = "bold",
  #   cat.col = c("red", "orange", "blue"),
  #   cat.cex = 0.8,
  #   cat.fontfamily = "serif"
  # )
  
  
  #return(re.out)
}

#' Title
#'
#' @param dir.name
#' @param input.file.pattern
#'
#' @return
#' @export
#'
#' @examples
#'
#' dir.name="/media/H_driver/2016/Yang/BedFromPeakCall/"
#' input.file.pattern="*bam_mm_summits_from_PeakCall4Yang.bed"
#'
#' out.dir.name="/media/H_driver/2016/Yang/Results/"
#' output.test.file= "SMC1A_RAD21_overlap_test"
#'
#'
#' dir.name.2="/media/H_driver/2016/Yang/BedFromPeakCall_mm_100e_07/"
#' input.file.pattern.2="*_bam_mm_1_from_PeakCall4Yang.bed"
#'
#' out.dir.name.2="/media/H_driver/2016/Yang/Results/"
#' output.test.file.2= "SMC1A_RAD21_mm_100e-07_peaks_overlap_test"
#'
#' ParserBedFile(dir.name,input.file.pattern,out.dir.name)
#'
#' re.mm.100e.07<-ParserBedFile(dir.name.2,input.file.pattern.2,out.dir.name)
#'
#'
#' save.image(file=paste0(out.dir.name,"re_save_2.RData"))
#' savehistory(file=paste0(out.dir.name,"re_save_2.Rhistory"))
#'
ParserBedFile<-function(dir.name,input.file.pattern,out.dir.name){
  
  file.name=paste0(dir.name,dir(dir.name,recursive = TRUE,pattern=input.file.pattern))
  file.name.2<-as.list(file.name)
  
  names(file.name.2)=sapply(strsplit(file.name,split="\\/"),"[[",7)
  
  print(file.name.2)
  
  file.name.3<-file.name.2[-6]
  
  re.out<-lapply(file.name.3,function(u){
    re=readBed(u)
    #colnames(re)=c("Count","GeneName")
    re
  })
  
  sample.name<-sapply(strsplit(sapply(strsplit(names(re.out),split="_i"),"[[",1),split="2470IUPUI_"),"[[",2)
  
  A=re.out[[2]] #SMC1
  B=re.out[[4]] #Rad21
  C=re.out[[1]] #ASXL1
  
  A.name=sample.name[2] #SMC1 name
  B.name=sample.name[4] #Rad21 name
  C.name=sample.name[1] #ASXL1 name
  
  #A <- createRandomRegions(nregions=50, length.mean=5000000, length.sd=3000000)
  #B <- c(A[1:25], createRandomRegions(nregions=25, length.mean=500000, length.sd=30000))
  #class(A)
  #class(B)
  
  
  TestOverlap <- function(A, B, A.name,B.name,out.dir.name) {
    NumOverLap<-numOverlaps(A, B, count.once=TRUE)
    pt <- overlapPermTest(A=A, B=B, ntimes=50)
    png(filename=paste0(out.dir.name,A.name,B.name,"_overlap_test.png"))
    plot(pt)
    dev.off()
    
    lz <- localZScore(pt=pt, A=A, B=B)
    png(filename=paste0(out.dir.name,A.name,B.name,"_localZScore.png"))
    plot(lz)
    dev.off()
  }
  
  TestOverlap(A,B,A.name,B.name,out.dir.name)
  TestOverlap(A,C,A.name,C.name,out.dir.name)
  TestOverlap(B,C,B.name,C.name,out.dir.name)
  
  # temp.name<-strsplit(names(file.name.2),split="\\.")
  
  #a = readBed(system.file("extdata", "examples/combined_regions.bed",
  #                        package="LOLA"))
  
  return(re.out)
}

#' ParserBedFile4Peng
#'
#' @param dir.name
#' @param input.file.pattern
#'
#' @return
#' @export
#'
#' @examples
#'
#' dir.name="/media/H_driver/2016/Yang/BedFromPeng/"
#' input.file.pattern="*.bed"
#'
#' out.dir.name="/media/H_driver/2016/Yang/Results/"

#'
#' re.from.bed.peng<-ParserBedFile4Peng(dir.name,input.file.pattern,out.dir.name)
#'
#'
#' save.image(file=paste0(out.dir.name,"re_save_peng.RData"))
#' savehistory(file=paste0(out.dir.name,"re_save_peng.Rhistory"))
#'
ParserBedFile4Peng<-function(dir.name,input.file.pattern,out.dir.name){
  
  file.name=paste0(dir.name,dir(dir.name,recursive = TRUE,pattern=input.file.pattern))
  file.name.2<-as.list(file.name)
  
  names(file.name.2)=sapply(strsplit(file.name,split="\\/"),"[[",7)
  
  print(file.name.2)
  
  #file.name.3<-file.name.2[-6]
  
  re.out<-lapply(file.name.2,function(u){
    re=readBed(u)
    #colnames(re)=c("Count","GeneName")
    re
  })
  
  sample.name<-sapply(strsplit(names(re.out),split="_peaks_"),"[[",1)
  
  A=re.out[[1]] #10_2470IUPUI_WT_BM_ASXL1
  B=re.out[[2]] #10_WT_BM_ASXL1
  C=re.out[[3]] #11_2470IUPUI_WT_BM_SMC1
  D=re.out[[5]] #13_2470IUPUI_WT_BM_Rad21
  
  A.name=sample.name[1] #10_2470IUPUI_WT_BM_ASXL1
  B.name=sample.name[2] #10_WT_BM_ASXL1
  C.name=sample.name[3] #11_2470IUPUI_WT_BM_SMC1
  D.name=sample.name[5] #13_2470IUPUI_WT_BM_Rad21
  
  #A <- createRandomRegions(nregions=50, length.mean=5000000, length.sd=3000000)
  #B <- c(A[1:25], createRandomRegions(nregions=25, length.mean=500000, length.sd=30000))
  #class(A)
  #class(B)
  
  
  TestOverlap <- function(A, B, A.name,B.name,out.dir.name,n_permutations) {
    NumOverLap<-numOverlaps(A, B, count.once=TRUE)
    pt <- overlapPermTest(A=A, B=B, ntimes=n_permutations)
    png(filename=paste0(out.dir.name,A.name,B.name,n_permutations,"_overlap_test.png"))
    plot(pt)
    dev.off()
    
    lz <- localZScore(pt=pt, A=A, B=B)
    png(filename=paste0(out.dir.name,A.name,B.name,n_permutations,"_localZScore.png"))
    plot(lz)
    dev.off()
  }
  
  TestOverlap(A,B,A.name,B.name,out.dir.name,100)
  TestOverlap(A,C,A.name,C.name,out.dir.name,100)
  TestOverlap(A,D,A.name,D.name,out.dir.name,100)
  
  TestOverlap(B,C,B.name,C.name,out.dir.name,100)
  TestOverlap(B,D,B.name,D.name,out.dir.name,100)
  TestOverlap(C,D,C.name,D.name,out.dir.name,100)
  
  
  # temp.name<-strsplit(names(file.name.2),split="\\.")
  
  #a = readBed(system.file("extdata", "examples/combined_regions.bed",
  #                        package="LOLA"))
  
  return(re.out)
}

#' ParserBedFile4PengVenn
#'
#' @param dir.name
#' @param input.file.pattern
#'
#' @return
#' @export
#'
#' @examples
#'
#' dir.name="/media/H_driver/2016/Yang/BedFromPeng/"
#' input.file.pattern="*combine_chr_st_end.bed"
#'
#' input.file.pattern="*combine_chr_st.bed"
#'
#'
#'
#' out.dir.name="/media/H_driver/2016/Yang/Results/"
#'
#' re.from.bed.peng.4.venn<-ParserBedFile4PengVenn(dir.name,input.file.pattern,out.dir.name)
#'
#'
#' save.image(file=paste0(out.dir.name,"re_save_peng.RData"))
#' savehistory(file=paste0(out.dir.name,"re_save_peng.Rhistory"))
#'
ParserBedFile4PengVenn<-function(dir.name,input.file.pattern,out.dir.name){
  
  file.name=paste0(dir.name,dir(dir.name,recursive = TRUE,pattern=input.file.pattern))
  file.name.2<-as.list(file.name)
  
  names(file.name.2)=sapply(strsplit(file.name,split="\\/"),"[[",7)
  
  print(file.name.2)
  
  #file.name.3<-file.name.2[-6]
  
  re.out<-lapply(file.name.2,function(u){
    re=read.table(u,header=FALSE)
    re<-as.character(re[,1])
    #colnames(re)=c("Count","GeneName")
    re
  })
  
  sample.name<-sapply(strsplit(names(re.out),split="_peaks_"),"[[",1)
  
  names(re.out)=sample.name
  
  A=re.out[[1]] #10_2470IUPUI_WT_BM_ASXL1
  B=re.out[[2]] #10_WT_BM_ASXL1
  C=re.out[[3]] #11_2470IUPUI_WT_BM_SMC1
  D=re.out[[5]] #13_2470IUPUI_WT_BM_Rad21
  
  A.name=sample.name[1] #10_2470IUPUI_WT_BM_ASXL1
  B.name=sample.name[2] #10_WT_BM_ASXL1
  C.name=sample.name[3] #11_2470IUPUI_WT_BM_SMC1
  D.name=sample.name[5] #13_2470IUPUI_WT_BM_Rad21
  
  #A <- createRandomRegions(nregions=50, length.mean=5000000, length.sd=3000000)
  #B <- c(A[1:25], createRandomRegions(nregions=25, length.mean=500000, length.sd=30000))
  #class(A)
  #class(B)
  
  
  # TestOverlap <- function(A, B, A.name,B.name,out.dir.name,n_permutations) {
  #   NumOverLap<-numOverlaps(A, B, count.once=TRUE)
  #   pt <- overlapPermTest(A=A, B=B, ntimes=n_permutations)
  #   png(filename=paste0(out.dir.name,A.name,B.name,n_permutations,"_overlap_test.png"))
  #   plot(pt)
  #   dev.off()
  #
  #   lz <- localZScore(pt=pt, A=A, B=B)
  #   png(filename=paste0(out.dir.name,A.name,B.name,n_permutations,"_localZScore.png"))
  #   plot(lz)
  #   dev.off()
  # }
  #
  # TestOverlap(A,B,A.name,B.name,out.dir.name,100)
  # TestOverlap(A,C,A.name,C.name,out.dir.name,100)
  # TestOverlap(A,D,A.name,D.name,out.dir.name,100)
  #
  # TestOverlap(B,C,B.name,C.name,out.dir.name,100)
  # TestOverlap(B,D,B.name,D.name,out.dir.name,100)
  # TestOverlap(C,D,C.name,D.name,out.dir.name,100)
  #
  
  # temp.name<-strsplit(names(file.name.2),split="\\.")
  
  #a = readBed(system.file("extdata", "examples/combined_regions.bed",
  #                        package="LOLA"))
  
  
  venn.plot <- venn.diagram(
    x = re.out[c(1,2)],
    filename = paste0(out.dir.name,names(re.out)[1],names(re.out)[2],"_overlap_venn.tiff"),
    height = 3000,
    width = 3500,
    resolution = 1000,
    col = "black",
    lty = "dotted",
    lwd = 1,
    fill = c("red","blue"),
    alpha = 0.50,
    label.col = c(rep("white",3)),
    cex = 0.5,
    fontfamily = "serif",
    fontface = "bold",
    cat.col = c("red","blue"),
    cat.cex = 0.5,
    cat.pos = 0.5,
    cat.dist = 0.05,
    cat.fontfamily = "serif"
  )
  
  venn.plot <- venn.diagram(
    x = re.out[c(1,3,5)],
    filename =paste0(out.dir.name,names(re.out)[1],names(re.out)[3],names(re.out)[5],"_overlap_venn.tiff"),
    col = "black",
    lty = "dotted",
    lwd = 2,
    fill = c("red", "orange", "blue"),
    alpha = 0.50,
    label.col = c(rep("white",7)),
    cex = 1,
    fontfamily = "serif",
    fontface = "bold",
    cat.col = c("red", "orange", "blue"),
    cat.cex = 0.8,
    cat.fontfamily = "serif"
  )
  
  
  return(re.out)
}

#' ParserBedFile4PengXls
#'
#' @param dir.name
#' @param input.file.pattern
#'
#' @return
#' @export
#'
#' @examples
#'
#' dir.name="/media/H_driver/2016/Yang/MACS/MACS/"
#' input.file.pattern="*.xls"
#'
#' out.dir.name="/media/H_driver/2016/Yang/Results/"

#'
#' re.from.bed.peng.Xls<-ParserBedFile4PengXls(dir.name,input.file.pattern,out.dir.name)
#'
#' re.from.bed.peng.Xls.2<-ParserBedFile4PengXls(dir.name,input.file.pattern,out.dir.name)
#'
#' re.from.bed.peng.Xls.3<-ParserBedFile4PengXls(dir.name,input.file.pattern,out.dir.name)
#'
#' save.image(file=paste0(out.dir.name,"re_save_peng.RData"))
#' savehistory(file=paste0(out.dir.name,"re_save_peng.Rhistory"))
#'
ParserBedFile4PengXls<-function(dir.name,input.file.pattern,out.dir.name){
  
  dir.name=reformatPath(dir.name)
  out.dir.name=reformatPath(out.dir.name)
  
  file.name=paste0(dir.name,dir(dir.name,recursive = TRUE,pattern=input.file.pattern))
  file.name.2<-as.list(file.name)
  
  names(file.name.2)=sapply(strsplit(file.name,split="\\/"),"[[",8)
  
  print(file.name.2)
  
  re.out<-lapply(names(file.name.2),function(u){
    
    
    #print(names(u))
    uu=file.name.2[[u]]
    
    re=read.table(uu,header=TRUE)
    colnames(re)[c(7,9)]=c("-10*LOG10(pvalue)","FDR(%)")
    temp<-re
    temp$start<-temp$start-1
    temp$end<-temp$end-1
    temp$summit<-temp$summit-1
    temp
    
    summitPeak<-temp$start+temp$summit
    
    temp2<-temp
    
    temp2$start<-summitPeak-49
    temp2$end<-summitPeak+50
    temp2$summit<-summitPeak
    temp2$length<-temp2$end-temp2$start+1
    temp2
    rownames(temp2)<-paste0("MACS_peak_",rownames(temp2))
    temp3<-toGRanges(temp2)
    names(temp3)<-rownames(temp2)
    
    
    dd.GRCm39.mm10<-toGRanges(EnsDb.Mmusculus.v75)
    genome(temp3)<-genome(dd.GRCm39.mm10)
    
    seqlevels(dd.GRCm39.mm10,force=TRUE) <- seqlevels(temp3)
    
    seqinfo(temp3)<-seqinfo(dd.GRCm39.mm10)
    
    temp3.trimmed<-trim(temp3, use.names=TRUE)
    
    genome.mm10<-getBSgenome("BSgenome.Mmusculus.UCSC.mm10")
    
    genome(temp3.trimmed)<-"mm10"
    
    seq.temp3<-getAllPeakSequence(temp3.trimmed,genome=genome.mm10)
    
    write2FASTA(seq.temp3, paste0(out.dir.name,u,".fa"))
    
    re2<-list(originalPeak=re,down1Peak=temp,aroundSummit100Peak=temp2,GR=temp3)
    re2
  })
  
  sample.name<-sapply(strsplit(names(file.name.2),split="_peaks_"),"[[",1)
  
  names(re.out)<-sample.name
  
  
  return(re.out)
}

#' ParserReadFiles
#'
#' @param input.file.dir 
#' @param output.file.dir 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' input.file.dir="~/TestBam/"
#' out.file.dir="~/"
#' testR<-ParserReadFiles(input.file.dir) 
#' 
#' 
ParserReadFiles <- function(input.file.dir,input.file.type) {
  
  #For input
  dir.name=input.file.dir
  dir.name=reformatPath(dir.name)
  
  file.name=file.path(dir.name,dir(dir.name,recursive = TRUE))
  
  file.name.2<-as.list(file.name)
  
  file.name.3<-lapply(1:length(file.name.2),function(u,input.file.type,file.name.2){
    
    tmp=file.name.2
    
    x=tmp[[u]]
    path_name=dirname(x)
    file_name=basename(x)
    
    n <- length(grep(input.file.type,file_name))
    
    #if(file_ext(file_name)==input.file.type)
    if(n == 1) {
      re<-file.path(path_name,file_name)
      #names(re)[u]<-file_name
    }else
    {re<-NULL}
    re
  },input.file.type,file.name.2)
  
  file.name.4<-file.name.3[lapply(file.name.3, length) > 0]
  
  names(file.name.4)=unlist(lapply(1:length(file.name.4),function(u,file.name.4){
    tmp=file.name.4
    x=tmp[[u]]
    path_name=dirname(x)
    file_name=basename(x)
    file_name
  },file.name.4))
  
  #print(file.name.4)
  
  #For output
  #output.dir.name=reformatPath(output.file.dir)
  #print(output.dir.name)
  
  #temp=Sys.time()
  #temp1=gsub(":","-",temp)
  #temp2=gsub(" ","-",temp1)
  #temp3=file.path(output.dir.name,"ReadBam_at_",temp2)
  #temp3=output.dir.name
  #dir.create(temp3)
  
  re2<-list(input=file.name.4)
  
  return(re2)
}

#' ChipSeq:::getGeneBedName("/Volumes/Bioinformatics$/2017/DannyNewData/NewRe2Danny","peaks.csv","/Volumes/Bioinformatics$/2017/DannyNewData/NewRe2Danny")
#' 
getGeneBedName <- function(input.file.dir,file.pattern,output.file.dir) {
  
  re<-ParserReadFiles(input.file.dir,file.pattern)
  
  file.name.2<-re$input
  
  re.out<-file.name.2
  
  if (!dir.exists(output.file.dir))
  {
    dir.create(output.file.dir, recursive = TRUE)
  }
  
  temp3 =output.file.dir
  
  cmd.1 <- lapply(1:length(re.out),function(u,re.out,temp3){
    
    x <- read.csv(re.out[[u]],header = TRUE)
    
    xx <- x[grep("Promoter",x$annotation),]
    
    chr <- which(colnames(xx) == "seqnames")
    gs  <- which(colnames(xx) == "geneStart")
    ge  <- which(colnames(xx) == "geneEnd")
    
    file_name = file_path_sans_ext(basename(re.out[[u]]))
    
    write.table(xx[,c(chr,gs,ge)],file=file.path(temp3,paste0(file_name,"_gene.bed")),quote = F,col.names = F,row.names = F,sep="\t")
    
    write.table(unique(xx$SYMBOL),file=file.path(temp3,paste0(file_name,"_gene_name.txt")),quote = F,col.names = F,row.names = F)
    
    ps  <- which(colnames(xx) == "start")
    pe  <- which(colnames(xx) == "end")
    
    write.table(xx[,c(chr,ps,pe)],file=file.path(temp3,paste0(file_name,"_peaks.bed")),quote = F,col.names = F,row.names = F,sep="\t")
    
  },re.out,temp3)
  
}

#' PeakCallAndAnnotation
#'
#' Call peak for bam files using macs14 and perform peak annotation using ChIPpeakAnno
#' 
#' @param input.file.dir 
#' @param input.file.pattern 
#' @param index.file 
#' @param output.file.dir 
#' @param genome 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' input.file.dir="/projects/scratch/bbc/Project/Danny_chip/Alignment/BWA/"
#' input.file.pattern="*.bam"
#' index.file=9
#' output.file.dir="/scratch/projects/bbc/aiminy_project/
#' genome="Hs"
#' 
#' PeakCallAndAnnotation(input.file.dir,output.file.dir,genome)
#' 
PeakCallAndAnnotation <- function(input.file.dir,output.file.dir,genome) {
  
  
  re<-ParserReadFiles(input.file.dir,"bam",output.file.dir)
  
  
  file.name.2<-re$input
  output.dir.name=re$output
  
  
  temp3=paste0(output.dir.name,"_PeakCall")
  
  dir.create(temp3)
  
  re.out<-file.name.2
  
  cmd9="macs14 -t "
  cmd10="-f BAM -g hs -n "
  cmd11=" -m 6,18 --bw=200 -p 0.00001"
  
  lapply(1:length(re.out),function(u,re.out,temp3){
    
    x=re.out[[u]]
    x_name=names(re.out)[u]
    
    cmd12=paste(cmd9,x,cmd10,paste0(temp3,"/",paste0(x_name,"_hs_1.00e-05_macs142")),cmd11,sep=" ")
    
    print(cmd12)
    system(cmd12, intern = TRUE, ignore.stderr = TRUE)
    
    #re=read.table(u,header=FALSE)
    #  re<-as.character(re[,1])
    #  #colnames(re)=c("Count","GeneName")
    #  re
  },re.out,temp3)
  
  #AnnotatePeak2(paste0(temp3,"/"),"*macs142_peaks.bed",7,paste0(output.dir.name,"PeakAnnotation_at_",temp2),genome="Hs")
  
  AnnotatePeak3(paste0(temp3,"/"),paste0(output.dir.name,"_PeakAnnotation"),
                genome="Hs")
  
  BamFileSortIndexVisualization2(re,genome)
  
}

BamFileSortIndexVisualization2 <- function(input,genome) {
  
  #library(ChIPpeakAnno)
  
  re<-input
  
  file.name.2<-re$input
  output.dir.name=re$output
  
  temp3=paste0(output.dir.name,"_visualization")
  
  dir.create(temp3)
  
  re.out<-file.name.2
  
  cmd1="samtools sort"
  
  lapply(1:length(re.out),function(u,re.out,temp3){
    
    x=re.out[[u]]
    x_name=names(re.out)[u]
    cmd2=paste0(cmd1," ",x," ",paste0(temp3,"/",x_name,"_sorted"))
    
    #print(cmd2)
    system(cmd2)
    
  },re.out,temp3)
  
  cmd3="samtools index"
  
  lapply(1:length(re.out),function(u,re.out,temp3){
    
    x=re.out[[u]]
    x_name=names(re.out)[u]
    cmd4=paste0(cmd3," ",paste0(temp3,"/",x_name,"_sorted.bam"))
    
    #print(cmd2)
    system(cmd4)
    
  },re.out,temp3)
  
  cmd5="ngs.plot.r -G hg19 -R tss -C"
  cmd6="-O"
  cmd7="-L 4000"
  #cmd3="-L 4000 -RR 1 -CD 1 -CO \\\"blue\\\""
  
  #ngs.plot.r -G hg19 -R tss -C $1 -O $2 -L 4000 -RR 1 -CD 1 -CO "blue"
  
  #file.name.3<-file.name.2[-6]
  
  lapply(1:length(re.out),function(u,re.out,temp3){
    
    x=re.out[[u]]
    x_name=names(re.out)[u]
    
    cmd8=paste(cmd5,paste0(temp3,"/",paste0(x_name,"_sorted.bam")),cmd6,paste0(temp3,"/",paste0(x_name,"_sorted")),cmd7,sep=" ")
    
    print(cmd8)
    system(cmd8, intern = TRUE, ignore.stderr = TRUE)
    
    #re=read.table(u,header=FALSE)
    #  re<-as.character(re[,1])
    #  #colnames(re)=c("Count","GeneName")
    #  re
  },re.out,temp3)
  
}

#' ProcessUsingCHIPpeakAnno
#'
#' @param input.file.dir
#' @param input.file.pattern
#' @param output.file.dir
#'
#' @return
#' @export
#'
#' @examples
#' input.file.dir="/media/H_driver/2016/Yang/MACS/MACS/"
#' input.file.pattern="*.bed"
#' output.file.dir="/media/H_driver/2016/Yang/MACS/MACS/"
#'
#' ProcessUsingCHIPpeakAnno(input.file.dir,input.file.pattern,output.file.dir)
#'
ProcessUsingCHIPpeakAnno <- function(input.file.dir,input.file.pattern,output.file.dir) {
  
  library(ChIPpeakAnno)
  
  dir.name=input.file.dir
  input.file.pattern=input.file.pattern
  
  dir.name=reformatPath(dir.name)
  
  output.dir.name=reformatPath(output.file.dir)
  dir.create(output.dir.name)
  
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
  
  library(TxDb.Hsapiens.UCSC.mm9.knownGene)
  
  aCR<-assignChromosomeRegion(overlaps, nucleotideLevel=FALSE,
                              precedence=c("Promoters", "immediateDownstream",
                                           "fiveUTRs", "threeUTRs",
                                           "Exons", "Introns"),
                              TxDb=TxDb.Mmusculus.UCSC.mm9.knownGene)
  
  
  
  
  barplot(aCR$percentage, las=3)
  
  pie1(aCR$percentage,las=3)
  
  dc<-annoGR(TxDb.Mmusculus.UCSC.mm9.knownGene)
  
  seqinfo(dc)
  seqlevels(dc)
  
  seqinfo(overlaps)
  seqlevels(overlaps)
  
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

# input.file.dir="/media/H_driver/2016/Danny/Danny_chip_PeakCall/"
# input.file.pattern="*macs142_peaks.bed"
# output.file.dir="/media/H_driver/2016/Danny/Danny_chip_PeakCall/"
#
# AnnotatePeak(input.file.dir,input.file.pattern,8,output.file.dir,genome="Mm")
# AnnotatePeak(input.file.dir,input.file.pattern,7,output.file.dir,genome="Hs")

AnnotatePeak <- function(input.file.dir,input.file.pattern,index.file,output.file.dir,genome) {
  
  library(ChIPpeakAnno)
  
  dir.name=input.file.dir
  input.file.pattern=input.file.pattern
  
  dir.name=reformatPath(dir.name)
  output.dir.name=reformatPath(output.file.dir)
  
  #print(output.dir.name)
  temp=Sys.time()
  temp1=gsub(":","-",Sys.time())
  temp2=gsub(" ","-",temp1)
  temp3=paste0(output.dir.name,"AnalysisResults_at_",temp2)
  
  dir.create(temp3)
  
  file.name=paste0(dir.name,dir(dir.name,recursive = TRUE,pattern=input.file.pattern))
  file.name.2<-as.list(file.name)
  
  names(file.name.2)=sapply(strsplit(file.name,split="\\/"),"[[",index.file)
  
  print(file.name.2)
  
  file.name.3<-file.name.2
  
  #sample.name<-sapply(strsplit(names(file.name.3),split="_peaks_"),"[[",1)
  
  #names(file.name.3)=sample.name
  
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
  
  if(genome=="Mm"){
    
    annoData <- toGRanges(EnsDb.Mmusculus.v75, feature="gene")
    
    ol <- findOverlapsOfPeaks(re.out[c(2,4,1)])
    
    overlaps<-ol$peaklist$`11_2470IUPUI_WT_BM_SMC1_peaks.bed///13_2470IUPUI_WT_BM_Rad21_peaks.bed///10_WT_BM_ASXL1_peaks.bed`
    
    binOverFeature(overlaps, annotationData=annoData,
                   radius=5000, nbins=20, FUN=length, errFun=0,
                   ylab="count",
                   main="Distribution of aggregated peak numbers around TSS")
    
    overlaps.trimmed<-trim(overlaps, use.names=TRUE)
    
    library(EnsDb.Mmusculus.v79)
    dd.GRCm39.mm10<-toGRanges(EnsDb.Mmusculus.v75)
    #seqinfo(dd.GRCm39.mm10)
    #seqlevels(dd.GRCm39.mm10)
    
    seqlevels(dd.GRCm39.mm10,force=TRUE) <- c("chr1","chr10","chr11","chr12","chr13",
                                              "chr14","chr15","chr16","chr17","chr18","chr19","chr2",
                                              "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY")
    
    seqinfo(overlaps)<-seqinfo(dd.GRCm39.mm10)
    #GRCm38/mm10
    #dd<-toGRanges(EnsDb.Mmusculus.v79)
    #seqinfo(dd)
    #library(ensembldb)
    #library(GenomeInfoDb)
    seqlevelsStyle(overlaps.trimmed) <- seqlevelsStyle(dd.GRCm39.mm10)
    overlaps.anno<-annoPeaks(overlaps.trimmed,dd.GRCm39.mm10)
    write.table(overlaps.anno,file=paste0(temp3,"/","annotation.txt"),row.names = FALSE,quote=FALSE,sep="\t")
  }else if(genome=="Hs"){
    
    library(EnsDb.Hsapiens.v75)
    #annoData<-toGRanges(EnsDb.Hsapiens.v75, feature="gene")
    
    dd.hs<-toGRanges(EnsDb.Hsapiens.v75)
    
    print(seqinfo(dd.hs))
    print(seqlevels(dd.hs))
    
    #print(seqlevels(dd.hs)[,1])
    
    #print(seqlevels(re.out[[1]])[,1])
    
    # seqlevels(dd.hs,force=TRUE) <- c("chr1","chr10","chr11","chr12","chr13",
    #                                           "chr14","chr15","chr16","chr17","chr18","chr19","chr2",
    #                                           "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY")
    
    #temp4=
    
    re.out.L<-lapply(1:length(re.out),function(u,re.out,dd.hs){
      
      x=re.out[[u]]
      x_name=names(re.out)[u]
      
      seqlevels(dd.hs,force=TRUE)<-seqinfo(x)@seqnames
      #print(seqinfo(re.out.trimmed))
      #print(seqlevels(re.out.trimmed))
      seqinfo(x)<-seqinfo(dd.hs)
      #GRCm38/mm10
      #dd<-toGRanges(EnsDb.Mmusculus.v79)
      #seqinfo(dd)
      #library(ensembldb)
      #library(GenomeInfoDb)
      seqlevelsStyle(x) <- seqlevelsStyle(dd.hs)
      re.out.trimmed<-trim(x, use.names=TRUE)
      overlaps.anno<-annoPeaks(re.out.trimmed,dd.hs)
      
      write.table(overlaps.anno,file=paste0(temp3,"/",x_name,"_annotation.txt"),row.names = FALSE,quote=FALSE,sep="\t")
    },re.out,dd.hs)
    
  }
  
}

ProcessUsingLOLA <- function() {
  #dbPath = system.file("extdata", "hg19", package="LOLA")
  
  dbPath = "/media/H_driver/2016/Yang/Results/data/groups/lab_bock/shared/resources/regions/LOLACore/mm10"
  
  regionDB = loadRegionDB(dbLocation=dbPath)
  
  #read triple overlap peaks
  regionSetA = readBed("/media/H_driver/2016/Yang/Results/WT_triple_overlap.bed")
  
  #read each peak
  regionSet.WT.ASXL1=readBed("/media/H_driver/2016/Yang/BedFromPeng/10_WT_BM_ASXL1_peaks_from_PeakCall4Yang.bed")
  regionSet.WT.SMC1=readBed("/media/H_driver/2016/Yang/BedFromPeng/11_2470IUPUI_WT_BM_SMC1_peaks_from_PeakCall4Yang.bed")
  regionSet.WT.Rad21=readBed("/media/H_driver/2016/Yang/BedFromPeng/13_2470IUPUI_WT_BM_Rad21_peaks_from_PeakCall4Yang.bed")
  
  
  db.841<-regionDB$regionGRL
  
  userSets = GRangesList(regionSet.WT.ASXL1,regionSet.WT.SMC1,regionSet.WT.Rad21)
  
  db.841.2<-c(db.841,userSets)
  
  userUniverse<-buildRestrictedUniverse(userSets)
  
  userUniverse.844<-buildRestrictedUniverse(db.841.2)
  
  re.SMC1A.RAD21.only<-makeVennDiagram(re.out[c(2,4)],NameOfPeaks=c("SMC1A", "RAD21"),totalTest=35000)
  
  ol.SMC1A.RAD21.only <- findOverlapsOfPeaks(re.out[c(2,4)])
  
  ol.SMC1A.RAD21.only.2<-ol.SMC1A.RAD21.only$peaklist$`11_2470IUPUI_WT_BM_SMC1_peaks.bed///13_2470IUPUI_WT_BM_Rad21_peaks.bed`
  
  checkUniverseAppropriateness(ol.SMC1A.RAD21.only.2,userUniverse.844, cores = 1,fast = FALSE)
  
  locResults2 = runLOLA(userSets,userUniverse,regionDB, cores=1)
  locResults = runLOLA(regionSetA,userUniverse,regionDB, cores=1)
  
  #SMC1A and RAD21 overlap
  
  
  
  locResults.SMC1A.RAD21 = runLOLA(ol.SMC1A.RAD21.only.2,userUniverse,regionDB, cores=4)
  
  head(locResults.SMC1A.RAD21)
  
  pValue<-as.data.frame(10^(-locResults.SMC1A.RAD21$pValueLog))
  
  temp<-as.data.frame(locResults.SMC1A.RAD21)
  
  locResults.SMC1A.RAD21.with.p<-cbind(temp[,1:3],pValue,temp[,4:23])
  
  colnames(locResults.SMC1A.RAD21.with.p)[4]="pValue"
  
  #It takes extensive computation resource, be carefully to run
  locResults.SMC1A.RAD21.based.large.universe = runLOLA(ol.SMC1A.RAD21.only.2,userUniverse.844,regionDB, cores=4)
  
  write.table(locResults.SMC1A.RAD21.with.p,
              file=paste0(out.dir.name,"SMC1A.RAD21.with.Other.Antibody.csv"),sep=",",
              quote = FALSE,row.names = FALSE,col.names = TRUE)
  
  write.table(locResults,file=paste0(out.dir.name,"Other.csv"),sep=",",quote = FALSE,row.names = FALSE,col.names = TRUE)
  
  listRegionSets(regionDB)
  extractEnrichmentOverlaps(locResults,regionSetA, regionDB)
  
  loadCaches("/media/H_driver/2016/Yang/Results/LOLACoreCaches_latest.tgz")
  
  save.image(file =paste0(out.dir.name,"re_search.RData"))
  
  # regionSetB = readBed("lola_vignette_data/setB_100.bed")
  # regionSetC = readBed("lola_vignette_data/setC_100.bed")
  # activeDHS = readBed("lola_vignette_data/activeDHS_universe.bed")
  # data("sample_universe", package="LOLA")
  #
  # dbPath = system.file("extdata", "hg19", package="LOLA")
  # regionDB = loadRegionDB(dbLocation=dbPath)
  # data("sample_universe", package="LOLA")
  # data("sample_input", package="LOLA")
  #
  # getRegionSet(regionDB, collections="ucsc_example", filenames="vistaEnhancers.bed")
  # getRegionSet(dbPath, collections="ucsc_example", filenames="vistaEnhancers.bed")
  #
  # res = runLOLA(userSets, userUniverse, regionDB, cores=1)
  # locResult = res[2,]
}

#' parser4diffbind
#' 
#' parser4diffbind
#'
#' @param input.sample.file
#' @param input.bed.dir
#' @param input.file.pattern
#'
#' @return
#' @export
#'
#' @examples
#'  
#' input.sample.file = "/Volumes/Bioinformatics$/2017/DannyNewData/SampleID_INFO_ChIP_new_Danny.csv"
#' input.bed.dir = "/Volumes/Bioinformatics$/2017/DannyNewData/PeakCall/"
#' input.file.pattern = "bed"
#' output.dir="BindDiff_4_7_2017-cutoff_10k_test"
#' output.dir="/Volumes/Bioinformatics$/2017/DannyNewData/testRun"
#' resmcf<-parser4diffbind(input.sample.file,input.bed.dir,input.file.pattern,output.dir)
#'
parser4diffbind<-function(input.sample.file,input.bed.dir,input.file.pattern,output.dir,
                          summit_peak=NULL){
  
  
  sample.info<-read.csv(input.sample.file,header = TRUE)
  sample.info.2 <- sample.info[,-dim(sample.info)[2]]
  
  out.dir.name = dirname(input.sample.file)
  
  re<-ParserReadFiles(input.bed.dir,input.file.pattern)
  
  re.bed<-re$input
  
  
  re.peaks.only.bed.2<-FL(re.bed,'peaks')
  
  if(!is.null(summit_peak)){
    re.summits.only.bed<-FL(re.bed,'summits')
  }
  
  #put peak files into database
  re.peaks.only.bed.3<-list_to_df(re.peaks.only.bed.2)
  
  
  print(re.peaks.only.bed.3)
  
  pos <- regexpr('R1', re.peaks.only.bed.3$ID)
  ID2 <- substr(re.peaks.only.bed.3$ID,1,pos+1)
  
  re.peaks.only.bed.3 <- cbind(ID2,re.peaks.only.bed.3)
  
  colnames(re.peaks.only.bed.3)=c("ID","ID2") 
  
  re.peaks.only.bed.4<-merge(sample.info.2,re.peaks.only.bed.3,by="ID")
  
  replicate<-apply(re.peaks.only.bed.4,1,function(x){
    ID=as.character(x[1])
    cell=as.character(x[2])
    TF=as.character(x[3])
    ID2=as.character(x[4])
    p1=nchar(ID)+2
    p2=regexpr(TF,ID2)
    p2=p2-2
    xx=substr(ID2,p1,p2)
    cell=substr(xx,1,nchar(xx)-2)
    rep=substr(xx,nchar(xx),nchar(xx))
    file=as.character(x[6])
    condition=paste(cell,TF,sep=":")
    z <- data.frame(t(c(ID2,cell,TF,condition,rep,file)))
    z
  })
  
  res <- do.call(rbind.data.frame,replicate)
  colnames(res) <- c("ID","Cell","Ab","Cell_Ab","Replicate","File")
  
  
  re.peaks.only.bed.5 <- res
  
  mcf<-NULL
  
  for(i in 1:dim(re.peaks.only.bed.5)[1]){
    
    #peaks.v=readBed(as.character(re.peaks.only.bed.5[i,6]))
    
    peaks.v <- readPeakFile(as.character(re.peaks.only.bed.5[i,6]), as = "GRanges")
    
    peak.caller.v="macs142"
    sampID.v=as.character(re.peaks.only.bed.5[i,1])
    
    sampID.v=substr(sampID.v,1,nchar(sampID.v)-5)
    
    tissue.v=as.character(re.peaks.only.bed.5[i,2])
    factor.v=as.character(re.peaks.only.bed.5[i,3])
    condition.v=as.character(re.peaks.only.bed.5[i,4])
    
    replicate.v=as.numeric(as.character(re.peaks.only.bed.5[i,5]))
    
    if(is.null(mcf)){
      mcf <- dba.peakset(NULL,peaks=peaks.v,peak.caller=peak.caller.v, sampID=sampID.v,tissue=tissue.v,
                         factor=factor.v,condition=condition.v,replicate=replicate.v)
    }else{
      mcf <- dba.peakset(mcf,peaks=peaks.v,peak.caller=peak.caller.v, sampID=sampID.v,tissue=tissue.v,
                         factor=factor.v,condition=condition.v,replicate=replicate.v)
    }
  }
  
  print(mcf)
  
  temp3=output.dir
  if(!dir.exists(temp3)){dir.create(temp3,recursive = TRUE)}
  
  save(mcf,file=file.path(temp3,"mcf.RData"))
  
  pdf(file.path(temp3,"corrmap.pdf"))
  dba.plotHeatmap(mcf,margin=23)
  dev.off()
  
  GetResultsFromDiffBind(mcf,"yes",temp3)
  
  return(mcf)
}

#' Title
#'
#' @param dir.name 
#'
#' @return
#' @export
#'
#' @examples
#' dir.name="/media/H_driver/2016/Yang/MACS/MACS/"
#' reformatPath(dir.name) 
reformatPath <- function(dir.name){
  
  CheckOPS=Sys.info()[['sysname']]
  
  if(CheckOPS=="Darwin"){
    
    if(length(grep('H_driver',dir.name))!=0){
      temp=unlist(strsplit(dir.name,split="\\/"))
      temp[2]="Volumes"
      temp[3]="Bioinformatics$"
      dir.name=paste0(paste0(temp,collapse = "/"),"/")
    }
  }
  
  return(dir.name)
}

#' useHOMER
#'
#' @return
#' @export
#'
#' @examples
#' res <- useHOMER("findMotifs.pl")
#' system(res)
#' 
useHOMER <- function(function.name){
  
  R_LIB=.libPaths()
  PATH1=Sys.getenv("PATH")
  
  homer_Lib=file.path(R_LIB,"ChipSeq/homer/bin")
  
  Sys.setenv(PATH=paste0(homer_Lib,":",PATH1)) 
  cmd1 <- Sys.which(function.name)[[1]]
  
  PATH2=Sys.getenv("PATH")
  
  cat(PATH1,"\n")
  cat(PATH2,"\n")
  
  return(cmd1)
}

# Filtering list based on certain pattern of the names of list elements

FL <- function(re.bed, patt) {
  
  re.peaks.only.bed<-lapply(1:length(re.bed),function(u,re.bed,patt){
    
    tmp=re.bed
    x=tmp[[u]]
    path_name=dirname(x)
    file_name=basename(x)
    
    #pos_summits=regexpr('summits',file_name)
    pos_peaks=regexpr(patt,file_name)
    
    if(pos_peaks>0){
      y<-x
    }else{
      y<-NULL  
    }
    
    y  
    
  },re.bed,patt)
  
  re.peaks.only.bed.2<-re.peaks.only.bed[lapply(re.peaks.only.bed,length) > 0]
  
  names(re.peaks.only.bed.2)=unlist(lapply(1:length(re.peaks.only.bed.2),function(u,re.peaks.only.bed.2){
    tmp=re.peaks.only.bed.2
    x=tmp[[u]]
    path_name=dirname(x)
    file_name=basename(x)
    file_name
  },re.peaks.only.bed.2))
  
  return(re.peaks.only.bed.2)
}

list_to_df <- function(list_for_df) {
  list_for_df <- as.list(list_for_df)
  
  nm <- names(list_for_df)
  
  if (is.null(nm)) 
    nm <- seq_along(list_for_df)
  
  ID<-sapply(1:length(nm),function(u,nm){
    x=nm[u]
    pos=regexpr("\\.",x)
    pos1=pos-1
    y<-substr(x,1,pos1)
    y
  },nm)
  
  df <- data.frame(ID=ID,name = nm, stringsAsFactors = FALSE)
  df$value <- unname(list_for_df)
  df
}

get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

#' plotBam2
#' 
#' @export
#' @example 
#' 
#' input.sample.file <- '/Volumes/Bioinformatics$/2017/DannyNewData/SampleID_INFO_ChIP_new_Danny.csv'
#' 
#' input.bam.file <- '/Users/axy148/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/ChipSeq/BamDanny/Bam_files.txt'
#' 
#' output.dir <- '/Users/axy148/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/ChipSeq/BamDanny/BamMerged'
#'  
#' re <- plotBam2(input.sample.file,input.bam.file,output.dir)
#' 
#' 
plotBam2 <- function(input.sample.file,input.bam.file,output.dir)
{
  
  if (!dir.exists(output.dir))
  {
    dir.create(output.dir, recursive = TRUE)
  }
  
  re <- GetSampleInfo(input.sample.file, input.bam.file)
  
  cellInfo <- re$re11
  
  cell.name <- unique(unlist(lapply(cellInfo$Type_Cell,function(u){substr(as.character(u),1,nchar(as.character(u))-2)})))
  
  cmd="samtools merge"
  
  cmd.3 <- lapply(cell.name,function(u,cellInfo,output.dir){
    
    cat(u,"\n")
    
    x <- cellInfo[which(substr(as.character(cellInfo$Type_Cell),1,nchar(as.character(cellInfo$Type_Cell))-2)==u),]
    x$Cell_TF <- gsub(" ","-",x$Cell_TF)
    
    tf <- as.character(x$Type_TF) 
    tf.unique <- unique(tf)
    
    y <- combn(1:length(tf),2,function(x)list(c(x[1],x[2])))
    
    yy <- lapply(tf.unique,function(u,tf,y){
      x <- which(tf %in% u)
      x
    },tf,y)
    
    yyy <- setdiff(y,yy)
    yyyy <- yyy[lapply(yyy, length) > 0]
    
    cmd.2 <- lapply(yyyy,function(u,x,output.dir){
      cmd.1 <- paste(cmd,file.path(output.dir,paste0(paste(x[c(u),]$Cell_TF,collapse ="-"),"_merged.bam")),paste(x[c(u),]$file.name,collapse = " "),collapse = " ")
    },x,output.dir)
    
    cmd.2
    
  },cellInfo,output.dir)
  
  cmd.4 <- unlist(cmd.3)
  
  mergeBamFilesUseJobArray(cmd.4)
  
}

mergeBamFilesUseJobArray <- function(cmd)
{
  
  index <- system('echo $LSB_JOBINDEX',intern = TRUE)
  total <- system('echo $LSB_JOBINDEX_END',intern = TRUE)
  
  cat(index,"\n\n")
  cat(total,"\n\n")
  
  u <- as.integer(index)
  
  system(cmd[u],intern = TRUE)
  cat(cmd,"\n\n")
  
}

#' R -e 'library(DoGs);library(ChipSeq);ChipSeq:::runPlotBam2("~/Danny/SampleID_INFO_ChIP_new_Danny.csv","~/Danny/Alignment/BWA/bam_files.txt","/scratch/projects/bbc/aiminy_project/DannyBamMerged")'

runPlotBam2 <-
  function(input.sample.file,
           input.bam.file,
           output.dir) {
    # This setting works
    Rfun1 <-
      'library(DoGs);library(ChipSeq);re <- ChipSeq:::plotBam2('
    
    input1 <- input.sample.file
    input2 <- input.bam.file
    output <- output.dir
    Rinput <-
      paste0('\\"',
             input1,
             '\\",',
             '\\"',
             input2,
             '\\",',
             '\\"',
             output,
             '\\"')
    Rfun2 <- ')'
    Rfun3 <- paste0(Rfun1, Rinput, Rfun2)
    
    n.job <- 12
    
    mergeBam <-
      DoGs:::createBsubJobArrayRfun(Rfun3, paste0("mergeBam[1-", n.job, "]"), NULL)
    
    system(mergeBam)
    
  }

#' bsub -P bbc -J "bamPlot" -o %J.bamPlot.log -e %J.bamPlot.err -W 72:00 -n 32 -q parallel -R 'rusage[mem= 16000 ] span[ptile= 16 ]' -u aimin.yan@med.miami.edu R -e 'library(ChipSeq);re <- ChipSeq:::plotMergedBamWgene("/scratch/projects/bbc/aiminy_project/DannyNewNgsPlot4Merged","*_sorted.bam","~/all_common_gene_unique.txt","NgsPlot4Merged")'
#'
plotMergedBamWgene <-
  function(input.bam.file.dir,
           file.pattern,
           input.gene.list,
           output.dir,
           ngs.para = c("hg19", 4000, 1, 1, "total"),
           add.input = NULL)
  {
    cellInfo <-
      list.files(
        input.bam.file.dir,
        pattern = file.pattern,
        all.files = TRUE,
        full.names = TRUE,
        recursive = TRUE,
        include.dirs = TRUE
      )
    
    
    if (!dir.exists(output.dir))
    {
      dir.create(output.dir, recursive = TRUE)
    }
    
    temp3 = output.dir
    
    cellInfo.run <- lapply(1:length(cellInfo), function(u, cellInfo,
                                                        temp3)
    {
      xx <- cellInfo[u]
      xxx <- file_path_sans_ext(basename(xx))
      cmd12 = paste(xx, input.gene.list, xxx, sep = "\t")
      cat(cmd12, file = file.path(temp3, paste0(xxx, "_config_cJunAndp27.txt")), sep = "\t")
    }, cellInfo, temp3)
    
    # dir.name=temp3 dir.name=reformatPath(dir.name)
    file.name <-
      list.files(
        temp3,
        pattern = "*_config_cJunAndp27.txt",
        all.files = TRUE,
        full.names = TRUE,
        recursive = TRUE,
        include.dirs = TRUE
      )
    
    file.name.2 <- as.list(file.name)
    
    zzz <- unlist(file.name.2)
    
    lapply(1:length(zzz), function(u, zzz)
    {
      dir.name = dirname(zzz[u][[1]])
      file_name = file_path_sans_ext(basename(zzz[u][[1]]))
      
      cmd = paste(
        "ngs.plot.r -G hg19 -R tss -C",
        zzz[u][[1]],
        "-O",
        file.path(dir.name,
                  paste0(file_name, ".tss")),
        "-T",
        file_name,
        "-L 4000 -RR 1 -CD 1 -GO total",
        sep = " "
      )
      
      cmd2 = cmd
      cat(cmd2, "\n")
      cat("\n")
      system(cmd2)
    }, zzz)
    
    return(re)
    
  }

#mcf.with.zhao <- ChipSeq:::addMoreBed2mcf("/Volumes/Bioinformatics$/2017/DannyNewData/SampleID_INFO_ChIP_new_Danny_zhao.csv","~/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/ChipSeq/BamDanny","/Volumes/Bioinformatics$/2017/DannyNewData/PeakCall","bed$",mcf=NULL,"/Volumes/Bioinformatics$/2017/DannyNewData/testRun5")

addMoreBed2mcf<-function(input.sample.file,input.bam.dir,input.bed.dir,input.file.pattern,mcf=NULL,output.dir){
  
  bam.file <- prepareBamFile(input.bam.dir)
  
  sample.info<-read.csv(input.sample.file,header = TRUE)
  sample.info.2 <- sample.info[,-dim(sample.info)[2]]
  
  #print(sample.info.2)
  
  #out.dir.name = dirname(input.sample.file)
  
  re<-ParserReadFiles(input.bed.dir,input.file.pattern)
  
  re.bed<-re$input
  
  #print(re.bed)
  
  re.peaks.only.bed.2<-Filter(function(x) !any(grepl("summits", x)), re.bed)
  
  # if(!is.null(summit_peak)){
  #   re.summits.only.bed<-FL(re.bed,'summits')
  # }
  
  #put peak files into database
  re.peaks.only.bed.3<-list_to_df(re.peaks.only.bed.2)
  
  print(re.peaks.only.bed.3)
  
  ll <- re.peaks.only.bed.3$ID
  
  ID2 <- unlist(lapply(ll,function(u) {
    
    pos <- regexpr('R1',u)
    if(pos>0){
      id <- substr(u,1,pos+1)
    }else{
      id <- u
    }
    id
  }))
  
  #   if(u>0){
  #   ID2 <- substr(ll,1,pos+1)
  #   }else{
  #     ID2 <- ll
  #   }
  
  # pos <- regexpr('R1', re.peaks.only.bed.3$ID)
  # 
  # print(pos)
  # 
  # ll <- re.peaks.only.bed.3$ID
  # 
  # ID2 <- lapply(pos,function(u,ll) {
  #  
  #   if(u>0){
  #   ID2 <- substr(ll,1,pos+1)
  #   }else{
  #     ID2 <- ll
  #   }
  #   
  # },ll)
  
  re.peaks.only.bed.3 <- cbind(ID2,re.peaks.only.bed.3)
  
  #print(re.peaks.only.bed.3)
  
  colnames(re.peaks.only.bed.3)=c("ID","ID2")
  #
  re.peaks.only.bed.4<-merge(sample.info.2,re.peaks.only.bed.3,by="ID")
  #
  
  cat("4\n")
  print(re.peaks.only.bed.4)
  cat("\n")
  
  re.peaks.only.bed.with.bam<-merge(bam.file,re.peaks.only.bed.4,by="ID")
  print(re.peaks.only.bed.with.bam)
  
  colnames(re.peaks.only.bed.with.bam)[6]="SampleID"
  
  replicate<-apply(re.peaks.only.bed.4,1,function(x){
    
    # pos <- regexpr('R1', re.peaks.only.bed.3$ID)
    
    ID=as.character(x[1])
    
    pos <- regexpr('R1',ID)
    if(pos>0){
      cell=as.character(x[2])
      TF=as.character(x[3])
      ID2=as.character(x[4])
      p1=nchar(ID)+2
      p2=regexpr(TF,ID2)
      p2=p2-2
      xx=substr(ID2,p1,p2)
      cell=substr(xx,1,nchar(xx)-2)
      rep=substr(xx,nchar(xx),nchar(xx))
      file=as.character(x[6])
      condition=paste(cell,TF,sep=":")
    }else
    {
      cell=as.character(x[2])
      TF=as.character(x[3])
      ID2=as.character(x[4])
      #p1=nchar(ID)+2
      #p2=regexpr(TF,ID2)
      #p2=p2-2
      #xx=substr(ID2,p1,p2)
      #cell=substr(xx,1,nchar(xx)-2)
      rep=1
      file=as.character(x[6])
      condition=paste(cell,TF,sep=":")
      
    }
    z <- data.frame(t(c(ID2,cell,TF,condition,rep,file)))
    z
  })
  
  res <- do.call(rbind.data.frame,replicate)
  colnames(res) <- c("ID","Cell","Ab","Cell_Ab","Replicate","File")
  
  colnames(res)[1] = "SampleID"
  
  sample.4.dba<-merge(res,re.peaks.only.bed.with.bam,by="SampleID")
  
  re.peaks.only.bed.5 <- res
  
  
  # replicate<-apply(re.peaks.only.bed.4,1,function(x){
  #   ID=as.character(x[1])
  #   cell=as.character(x[2])
  #   TF=as.character(x[3])
  #   ID2=as.character(x[4])
  #   p1=nchar(ID)+2
  #   p2=regexpr(TF,ID2)
  #   p2=p2-2
  #   xx=substr(ID2,p1,p2)
  #   cell=substr(xx,1,nchar(xx)-2)
  #   #rep=1
  #   file=as.character(x[6])
  #   condition=paste(cell,TF,sep=":")
  #   z <- data.frame(t(c(ID2,cell,TF,condition,rep,file)))
  #   z
  # })
  #
  # res <- do.call(rbind.data.frame,replicate)
  # colnames(res) <- c("ID","Cell","Ab","Cell_Ab","Replicate","File")
  #
  #
  # re.peaks.only.bed.5 <- res
  
  print(re.peaks.only.bed.5)
  
  mcf<-mcf
  
  ll <- list()
  
  for(i in 1:dim(re.peaks.only.bed.5)[1]){
    
    #peaks.v=readBed(as.character(re.peaks.only.bed.5[i,6]))
    
    peaks.v <- readPeakFile(as.character(re.peaks.only.bed.5[i,6]), as = "GRanges")
    peaks.v.sorted <- DiffBind:::pv.peaksort(as.data.frame(peaks.v))
    #peaks.v.sorted <- peaks.v
    ##print(peaks.v)
    ll[[i]] <- peaks.v.sorted
    peak.caller.v="macs142"
    sampID.v=as.character(re.peaks.only.bed.5[i,1])
    
    pos <- regexpr('R1',sampID.v)
    if(pos>0){
      sampID.v=substr(sampID.v,1,nchar(sampID.v)-5)
    }else{
      sampID.v = sampID.v
    }
    
    tissue.v=as.character(re.peaks.only.bed.5[i,2])
    factor.v=as.character(re.peaks.only.bed.5[i,3])
    condition.v=as.character(re.peaks.only.bed.5[i,4])
    
    replicate.v=as.numeric(as.character(re.peaks.only.bed.5[i,5]))
    
    if(is.null(mcf)){
      mcf <- dba.peakset(NULL,peaks=peaks.v.sorted,peak.caller=peak.caller.v, sampID=sampID.v,tissue=tissue.v,
                         factor=factor.v,condition=condition.v,replicate=replicate.v)
    }else{
      mcf <- dba.peakset(mcf,peaks=peaks.v.sorted,peak.caller=peak.caller.v, sampID=sampID.v,tissue=tissue.v,
                         factor=factor.v,condition=condition.v,replicate=replicate.v)
    }
  }
  
  print(mcf)
  
  temp3=output.dir
  if(!dir.exists(temp3)){dir.create(temp3,recursive = TRUE)}
  
  save(mcf,file=file.path(temp3,"mcf.RData"))
  
  pdf(file.path(temp3,"corrmap.pdf"))
  dba.plotHeatmap(mcf,margin=23)
  dev.off()
  
  # #GetResultsFromDiffBind(mcf,"yes",temp3)
  re <- list(mcf=mcf,ll = ll,re.peaks.only.bed.5 = re.peaks.only.bed.5,re.peaks.only.bed.with.bam = re.peaks.only.bed.with.bam, sample.4.dba =  sample.4.dba)
  return(re)
}

#ChipSeq:::GetResultsFromDiffBind2(temp,"/Volumes/Bioinformatics$/2017/DannyNewData/testRun2")

GetResultsFromDiffBind2<-function(mcf7,output.file.dir){
  
  temp<-mcf.with.zhao
  
  sampID.v<-colnames(temp$class) 
  sampID.v.2<-unlist(lapply(1:length(sampID.v),function(u,sampID.v){
    
    x=sampID.v[u]
    
    y=x
    y
  },sampID.v))
  
  
  colnames(temp$class)<-sampID.v.2
  temp$class[1,]<-sampID.v.2
  
  #temp<-dba(temp,mask =c(1,2,13,14,11,12,15,16))
  
  #Merge the replicates of each set 
  #if(Mergereplicates=="yes"){
  temp22<-dba.peakset(temp,consensus = -DBA_REPLICATE)
  #temp22<-dba(temp2,mask =c(9,16,17:23))
  print(temp22)
  ##need to set the flexiable number identified from data sets
  #    t <- temp2
  #    temp22<-dba(t,mask =which(t$class[which(rownames(t$class)=="Replicate"),]=="1-2"))
  #  }else{
  #    temp22<-temp
  #  }
  
  temp22<-dba(temp22,mask =17:28)
  print(temp22)
  
  po1 <- c(1,2,6,10)
  po2 <- c(1,2)
  
  temp3=file.path(output.file.dir,"Venn2")
  
  if(!dir.exists(temp3))
  {
    dir.create(temp3)
  }
  
  png(file.path(temp3,paste0(paste0(colnames(temp22$class)[po1],collapse = "-vs-"),".png")))
  dba.plotVenn(temp22,mask=po1,main=paste0(colnames(temp22$class)[po1],collapse = "-vs-"))
  dev.off()
  
  png(file.path(temp3,paste0(paste0(colnames(temp22$class)[po2],collapse = "-vs-"),".png")))
  dba.plotVenn(temp22,mask=po2,main=paste0(colnames(temp22$class)[po2],collapse = "-vs-"))
  dev.off()
  
  # Tissue1<-temp22$class[2,]
  # Tissue2<-unique(temp22$class[2,])
  # 
  # TF<-unique(temp22$class[3,])
  # TF.n<-length(TF)
  # 
  # temp3=file.path(output.file.dir,"Venn")
  # 
  # if(!dir.exists(temp3))
  # {
  #   dir.create(temp3)
  # }
  # 
  # for(i in 1:length(Tissue2)){
  #   po<-which(Tissue1 %in% Tissue2[i])
  #   
  #   print(po)
  #   
  #   if(length(po)==2)
  #   {
  #     png(file.path(temp3,paste0(paste0(colnames(temp22$class)[po],collapse = "-vs-"),".png")))
  #     dba.plotVenn(temp22,mask=po,main=Tissue2[i])
  #     dev.off()
  #     
  #   }else if(length(po)==4){
  #     po1<-po[c(1,3)]
  #     po2<-po[c(2,4)]
  #     
  #     png(file.path(temp3,paste0(paste0(colnames(temp22$class)[po1],collapse = "-vs-"),".png")))
  #     dba.plotVenn(temp22,mask=po1,main=Tissue2[i])
  #     dev.off()
  #     
  #     png(file.path(temp3,paste0(paste0(colnames(temp22$class)[po2],collapse = "-vs-"),".png")))
  #     dba.plotVenn(temp22,mask=po2,main=Tissue2[i])
  #     dev.off()
  #     
  #   }else
  #   {
  #     cat(paste0("For ",Tissue2[i],": Only one peak profile for ",TF.n," TFs\n"))
  #   }
  # }
  # 
  # p.common<-lapply(1:length(Tissue2),function(u,Tissue1,Tissue2,temp22){
  #   
  #   po<-which(Tissue1 %in% Tissue2[u])
  #   
  #   if(length(po)==2)
  #   {
  #     common.peaks<-dba.overlap(temp22,mask=po)
  #     y<-common.peaks$inAll
  #   }else if(length(po)==4){
  #     po1<-po[c(1,3)]
  #     po2<-po[c(2,4)]
  #     
  #     common.peaks.1<-dba.overlap(temp22,mask=po1)
  #     y1<-common.peaks.1$inAll
  #     
  #     common.peaks.2<-dba.overlap(temp22,mask=po2)
  #     y2<-common.peaks.2$inAll
  #     
  #     y<-list(y1=y1,y2=y2)
  #     
  #     names(y)[1]<-paste0(colnames(temp22$class)[po1],collapse = "-vs-")
  #     names(y)[2]<-paste0(colnames(temp22$class)[po2],collapse = "-vs-")
  #     
  #   }else{
  #     y<-NULL}
  #   y
  # },Tissue1,Tissue2,temp22)
  # 
  # names(p.common)<-Tissue2
  # 
  # p.common<-unlist(p.common,recursive = F)
  # 
  # p.common<-p.common[lapply(p.common,length) > 0]
  # 
  # if(!dir.exists(file.path(output.file.dir,"common_peaks_bed")))
  # {
  #   dir.create(file.path(output.file.dir,"common_peaks_bed"))
  #   dir.create(file.path(output.file.dir,"common_peaks_bed","ucsc"))
  #   dir.create(file.path(output.file.dir,"common_peaks_bed","igv"))
  # }
  # 
  # #output common peaks to bed files
  # 
  # lapply(1:length(p.common),function(u,p.common,output.file.dir){
  #   x=p.common[[u]]
  #   
  #   x_name=names(p.common)[u]
  #   
  #   df <- data.frame(seqnames=seqnames(x),
  #                    #starts=start(x)-1,
  #                    starts=start(x),
  #                    ends=end(x),
  #                    names=c(rep(".", length(x))),
  #                    scores=elementMetadata(x)[,1],
  #                    strands=strand(x))
  #   
  #   #assign strand
  #   df.str <- data.frame(seqnames=seqnames(x),
  #                        #starts=start(x)-1,
  #                        starts=start(x),
  #                        ends=end(x),
  #                        names=c(rep(".", length(x))),
  #                        scores=elementMetadata(x)[,1],
  #                        strands=c(rep(".", length(x))))
  #   
  #   df.str.1<-df.str[-grep("random",df.str$seqnames),]
  #   
  #   df.str.2<-df.str.1
  #   
  #   df.str.3<-df.str.2[-grep("chrUn",df.str.2$seqnames),]
  #   
  #   write.table(df,file=file.path(output.file.dir,"common_peaks_bed",paste0(x_name,"_cp_with_header.bed")),
  #               col.names=TRUE,row.names = FALSE,quote=FALSE,sep="\t")
  #   
  #   write.table(df,file=file.path(output.file.dir,"common_peaks_bed","ucsc",paste0(x_name,"_4_ucsc.bed")),
  #               col.names=FALSE,row.names = FALSE,quote=FALSE,sep="\t")
  #   
  #   write.table(df,file=file.path(output.file.dir,"common_peaks_bed",paste0(x_name,"_common_peaks.bed")),
  #               col.names=FALSE,row.names = FALSE,quote=FALSE,sep="\t")
  #   
  #   write.table(df,file=file.path(output.file.dir,"common_peaks_bed","igv",paste0(x_name,"_4_igv.bed")),
  #               col.names=FALSE,row.names = FALSE,quote=FALSE,sep="\t")
  #   
  # },p.common,output.file.dir)
  # 
  # AnntationUsingChipSeeker(file.path(output.file.dir,"common_peaks_bed","igv"),"bed",file.path(output.file.dir,"common_peaks_bed")
  #                          ,txdb="hg19",DD=5000,distanceToTSS_cutoff=10000)
  # 
  # return(p.common)
}

useChIPpeakAnno <- function(mcf.with.zhao) {
  #findOverlapsOfPeaks(mcf.with.zhao$ll[c(17,11,18)])
  overlap <- findOverlapsOfPeaks(mcf.with.zhao$ll[17][[1]],mcf.with.zhao$ll[11][[1]],mcf.with.zhao$ll[18][[1]])
  
  bed.li<- overlap$overlappingPeaks[2][[1]][,2:4]
  write.table(bed.li,file="cJun_Li.bed",sep="\t",row.names = FALSE,
              col.names = FALSE,quote=FALSE)
  
  bed.zhao <- overlap$overlappingPeaks[2][[1]][,8:10]
  write.table(bed.zhao,file="cJun_zhao.bed",sep="\t",row.names = FALSE,
              col.names = FALSE,quote=FALSE)
  
  s1 <- makeVennDiagram(mcf.with.zhao$ll[c(17,11,18)], NameOfPeaks=c("cJun_Li","cJun_1833","cJun_zhao"),
                        totalTest=20000,scaled=F, euler.d=F,fill = c("red","green","blue"),
                        alpha = 0.50,
                        label.col = c(rep("black",7)),
                        cex = 2,
                        fontfamily = "serif",
                        fontface = "bold",
                        cat.col = c("red","green","blue"),connectedPeaks = "keepAll")
  
  s2 <- makeVennDiagram(mcf.with.zhao$ll[c(17,3,18)], NameOfPeaks=c("cJun_Li","cJun_231","cJun_zhao"),
                        totalTest=20000,scaled=F, euler.d=F,fill = c("red","green","blue"),
                        alpha = 0.50,
                        label.col = c(rep("black",7)),
                        cex = 2,
                        fontfamily = "serif",
                        fontface = "bold",
                        cat.col = c("red","green","blue"),connectedPeaks = "keepAll")
}

# SampleID,Tissue,Factor,Condition,Treatment,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller
# BT4741,BT474,ER,Resistant,Full-Media,1,reads/Chr18_BT474_ER_1.bam,BT474c,reads/Chr18_BT474_input.bam,peaks/BT474_ER_1.bed.gz,bed
# BT4742,BT474,ER,Resistant,Full-Media,2,reads/Chr18_BT474_ER_2.bam,BT474c,reads/Chr18_BT474_input.bam,peaks/BT474_ER_2.bed.gz,bed
# MCF71,MCF7,ER,Responsive,Full-Media,1,reads/Chr18_MCF7_ER_1.bam,MCF7c,reads/Chr18_MCF7_input.bam,peaks/MCF7_ER_1.bed.gz,bed
# MCF72,MCF7,ER,Responsive,Full-Media,2,reads/Chr18_MCF7_ER_2.bam,MCF7c,reads/Chr18_MCF7_input.bam,peaks/MCF7_ER_2.bed.gz,bed
# MCF73,MCF7,ER,Responsive,Full-Media,3,reads/Chr18_MCF7_ER_3.bam,MCF7c,reads/Chr18_MCF7_input.bam,peaks/MCF7_ER_3.bed.gz,bed
# T47D1,T47D,ER,Responsive,Full-Media,1,reads/Chr18_T47D_ER_1.bam,T47Dc,reads/Chr18_T47D_input.bam,peaks/T47D_ER_1.bed.gz,bed
# T47D2,T47D,ER,Responsive,Full-Media,2,reads/Chr18_T47D_ER_2.bam,T47Dc,reads/Chr18_T47D_input.bam,peaks/T47D_ER_2.bed.gz,bed
# MCF7r1,MCF7,ER,Resistant,Full-Media,1,reads/Chr18_TAMR_ER_1.bam,TAMRc,reads/Chr18_TAMR_input.bam,peaks/TAMR_ER_1.bed.gz,bed
# MCF7r2,MCF7,ER,Resistant,Full-Media,2,reads/Chr18_TAMR_ER_2.bam,TAMRc,reads/Chr18_TAMR_input.bam,peaks/TAMR_ER_2.bed.gz,bed
# ZR751,ZR75,ER,Responsive,Full-Media,1,reads/Chr18_ZR75_ER_1.bam,ZR75c,reads/Chr18_ZR75_input.bam,peaks/ZR75_ER_1.bed.gz,bed
# ZR752,ZR75,ER,Responsive,Full-Media,2,reads/Chr18_ZR75_ER_2.bam,ZR75c,reads/Chr18_ZR75_input.bam,peaks/ZR75_ER_2.bed.gz,bed

#prepareBamFile("~/Dropbox (BBSR)/BWA")

prepareBamFile <- function(input.bam.dir,file.pattern){
  
  bam <- list.files(path = input.bam.dir, pattern=file.pattern, all.files = TRUE,
                    full.names = TRUE, recursive = TRUE,
                    ignore.case = FALSE, include.dirs = TRUE, no.. = FALSE)
  bam2 <- as.data.frame(cbind(basename(bam),bam))
  
  ID <- unlist(lapply(bam2[,1],function(u)
  {
    y <- u
    pos <- regexpr("\\.", y)
    pos <- pos - 1
    y <- substr(y, 1, pos)
    y
  }))
  
  bam3 <- cbind(ID,bam2)
  bam3
  
}

# tamoxifen <- prepareSampleSheet(mcf.with.zhao,"~/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/ChipSeq/BamDanny")

prepareSampleSheet <- function(mcf.with.zhao,output.file.dir){
  
  sampleSheet <- mcf.with.zhao$sample.4.dba[,c(1,2,3,5,9,13)]
  sampleSheet2 <- cbind(sampleSheet[,1:3],rep("Danny-Treatment",dim(sampleSheet)[1]),rep("Danny-condition",dim(sampleSheet)[1]),sampleSheet[,4:6],rep("bed",dim(sampleSheet)[1]))
  colnames(sampleSheet2)=c("SampleID","Tissue","Factor","Condition","Treatment","Replicate","bamReads","Peaks","PeakCaller")
  fwrite(sampleSheet2, file ="sampleSheet2.csv")
  tamoxifen <- dba(sampleSheet="sampleSheet2.csv")
  tamoxifen
  tamoxifen <- dba.count(tamoxifen, summits=250)
  
  tamoxifen <- dba.contrast(tamoxifen, group1 = c(1,2),group2=c(3,4),name1="cJun_231_DD", name2="cJun_231")
  tamoxifen <- dba.contrast(tamoxifen, group1 = c(1,2),group2=c(11,12),name1="cJun_231_DD", name2="cJun_1833")
  tamoxifen <- dba.contrast(tamoxifen, group1 = c(3,4),group2=c(11,12),name1="cJun_231", name2="cJun_1833")
  tamoxifen <- dba.contrast(tamoxifen, group1 = c(13,14),group2=c(5,6),name1="p27_231_DD", name2="p27_231")
  tamoxifen <- dba.contrast(tamoxifen, group1 = c(13,14),group2=c(15,16),name1="p27_231_DD", name2="p27_1833")
  tamoxifen <- dba.contrast(tamoxifen, group1 = c(5,6),group2=c(15,16),name1="p27_231", name2="p27_1833")
  tamoxifen <- dba.analyze(tamoxifen)
  
  save(tamoxifen,file=file.path(output.file.dir,"Danny_Dba.RData"))
  tamoxifen

  # dba.report(tamoxifen,contrast = 1)
  # dba.report(tamoxifen,contrast = 2)
  # dba.report(tamoxifen,contrast = 3)
  # dba.report(tamoxifen,contrast = 4)
  # dba.report(tamoxifen,contrast = 5)
  # dba.report(tamoxifen,contrast = 6)
  
}

#generateBedFromDba(tamoxifen,"~/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/ChipSeq/Bba2Bed")

generateBedFromDba <- function(tamoxifen,output.file.dir){
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  n <- length(tamoxifen$contrasts)
  
  lapply(1:n,function(u,tamoxifen){
    
    temp <- as.data.frame(dba.report(tamoxifen,contrast = u,th=1))
    
    file.name <- paste(colnames(as.data.frame(dba.report(tamoxifen,contrast = u,th=1)))[c(7,8)],collapse = "-")
    
    write.table(temp, row.names=F, col.names = F, quote =F,sep="\t", file=file.path(output.file.dir,paste0(file.name,".bed")))
    
  },tamoxifen)
  
}

#' dir.name="~/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/ChipSeq/Bba2Bed"
#' input.file.pattern=".bed"
#' out.dir.name="~/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/ChipSeq/AnnotationNew4DBA"
#' txdb="hg19"
#' DD=5000
#' 
#' AnntationUsingChipSeeker2(dir.name,input.file.pattern,out.dir.name,txdb=txdb,DD,distanceToTSS_cutoff=5000, AP=c("Promoter","Intron"))
#'
#' res.promoter <- AnntationUsingChipSeeker2(dir.name,input.file.pattern,out.dir.name,txdb=txdb,DD,distanceToTSS_cutoff=5000,AP=c("Promoter"))
#' 
#' AnntationUsingChipSeeker2(dir.name,input.file.pattern,out.dir.name,txdb=txdb,DD,distanceToTSS_cutoff=5000,AP=c("Intron"))
#' 
 
AnntationUsingChipSeeker2 <- function(dir.name,input.file.pattern,out.dir.name,txdb=c("hg19","hg38"),DD,distanceToTSS_cutoff=5000,assignGenomicAnnotation=TRUE,AP=c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic")) {
  
  re<-ParserReadFiles(dir.name,input.file.pattern)
  
  re.bed<-re$input
  
  re.peaks.only.bed.2 <- re.bed
  
  # if(length(dir(dir.name,pattern="peaks.bed"))!=0)
  # {
  # re.peaks.only.bed.2<-FL(re.bed,'peaks')
  # cat("peaks\n")
  # print(re.peaks.only.bed.2)
  # }
  # 
  # if(length(dir(dir.name,pattern="summits.bed"))!=0){
  # re.summits.only.bed<-FL(re.bed,'summits')
  # cat("summits\n")
  # print(re.summits.only.bed)
  # }
  
  txdb<-match.arg(txdb)
  
  switch (txdb,
          hg38 = {
            cat("Use hg38\n")
            txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
          },
          {
            cat("Use hg19\n") 
            txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
          }
  )
  
  APpath <- paste(AP,collapse = "_")
  
  temp3=file.path(out.dir.name,"Annotation",APpath)
  
  if(!dir.exists(temp3)){dir.create(temp3,recursive = TRUE)}
  
  d=DD
  
  peaks.anno.list <- lapply(1:length(re.peaks.only.bed.2),function(u,re.peaks.only.bed.2,d){
    
    peaks=readPeakFile(re.peaks.only.bed.2[[u]],as="data.frame")
    
    print(head(peaks))
    
    peakAnno <- annotatePeak(re.peaks.only.bed.2[[u]], tssRegion=c(-d, d),
                             TxDb=txdb,assignGenomicAnnotation=assignGenomicAnnotation,genomicAnnotationPriority=AP,annoDb="org.Hs.eg.db")
    
    #select.index <- which(peakAnno$distanceToTSS<20,000 && peakAnno$distanceToTSS >= -20,000)
    dropAnnoM <- function (csAnno, distanceToTSS_cutoff) 
    {
      idx <- which(abs(mcols(csAnno@anno)[["distanceToTSS"]]) < 
                     distanceToTSS_cutoff)
      csAnno@anno <- csAnno@anno[idx]
      csAnno@peakNum <- length(idx)
      if (csAnno@hasGenomicAnnotation) {
        csAnno@annoStat <- ChIPseeker:::getGenomicAnnoStat(csAnno@anno)
        csAnno@detailGenomicAnnotation = csAnno@detailGenomicAnnotation[idx, 
                                                                        ]
      }
      csAnno
    }
    
    peakAnno <- dropAnnoM(peakAnno,distanceToTSS_cutoff = distanceToTSS_cutoff)
    
    x_name=names(re.peaks.only.bed.2)[u]
    cat(x_name,"\n")
    png(file.path(temp3,paste0(x_name,"_",d,APpath,"_around_tss_annotation_pie.png")))
    plotAnnoPie(peakAnno)
    dev.off()
    
    peaks.anno=as.data.frame(peakAnno)
    
    print(head(peaks.anno))
    
    #print(paste0(peaks[,c(2,3)]))
    group1 <- strsplit(tools::file_path_sans_ext(x_name),"-")[[1]][1] 
    group2 <- strsplit(tools::file_path_sans_ext(x_name),"-")[[1]][2] 
    
    colnames(peaks.anno)[5:13]=c("starnd","width_DiffDind_based","strand","Conc", group1,group2,"Fold","p.value","FDR")
    
    print(colnames(peaks.anno))
    write.table(peaks.anno,file=file.path(temp3,paste0(x_name,"_",d,APpath,"_around_tss_annotation_4_only_mapped_peaks.xls")),
                row.names = FALSE,quote=FALSE,sep="\t")
    
    # unmapped.peaks<-peaks[-which(paste0(peaks[,2],"_",peaks[,3]) %in% paste0(peaks.anno[,2],"_",peaks.anno[,3])),]
    # 
    # cat(dim(peaks)[1]," ",dim(peaks.anno)[1]," ",dim(unmapped.peaks)[1],"\n")
    # 
    # 
    # if(dim(unmapped.peaks)[1]!=0){
    #   
    #   colnames(unmapped.peaks)=colnames(peaks.anno)[1:6]
    #   
    #   unmapped.peaks.3<-smartbind(peaks.anno,unmapped.peaks)
    #   
    #   unmapped.peaks.4<-unmapped.peaks.3[order(unmapped.peaks.3[,1],unmapped.peaks.3[,2]),]
    #   
    #   write.table(unmapped.peaks.4,file=file.path(temp3,paste0(x_name,"_",d,APpath,"_around_tss_annotation_4_all_peaks.xls")),row.names = FALSE,quote=FALSE,sep="\t")
    # }
    # 
    
    peaks.anno
    
  },re.peaks.only.bed.2,d)
  
}

#' out.dir.name="~/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/ChipSeq/CombineChipSeqWithRNASeq"
#' rna.seq.input.file="/Volumes/Bioinformatics$/2017/DannyNewData/AnnotationNew4DBA/table_S3_headerless.xlsx"
#' 
#' integrateChipSeqWithRNASeq(out.dir.name,rna.seq.input.file) 
#' 
integrateChipSeqWithRNASeq <- function(output.file.dir,rna.seq.input.file) {
 
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  gene.231dd.231.1833.1833shp27 <- read.xlsx(rna.seq.input.file,sheetIndex = 1,header = FALSE)
  
  colnames(gene.231dd.231.1833.1833shp27)=c("SYMBOL","FC_231DD_231","p_231DD_231","FC_1833_231","p_1833_231","FC_1833shp27_1833","p_1833shp27_1833")
  
  data.231DD.231 <- gene.231dd.231.1833.1833shp27[,c(1,2,3)]
  data.1833.231 <- gene.231dd.231.1833.1833shp27[,c(1,4,5)]
  data.1833shp27.1833 <- gene.231dd.231.1833.1833shp27[c(1,6,7)]
  
  cJun.231DD.231.ChipSeq.RNASeq <- merge(res.promoter[[2]],data.231DD.231,by="SYMBOL")
  p27.231DD.231.ChipSeq.RNASeq <- merge(res.promoter[[5]],data.231DD.231,by="SYMBOL")
  
  cJun.231.1833.ChipSeq.RNASeq <- merge(res.promoter[[3]],data.1833.231,by="SYMBOL")
  p27.231.1833.ChipSeq.RNASeq  <- merge(res.promoter[[6]],data.1833.231,by="SYMBOL")
  
  write.table(cJun.231DD.231.ChipSeq.RNASeq,file=file.path(output.file.dir,"cJun_231DD_231_DBA_DGE.xls"),
              row.names = FALSE,quote=FALSE,sep="\t")
  
  write.table(p27.231DD.231.ChipSeq.RNASeq,file=file.path(output.file.dir,"p27_231DD_231_DBA_DGE.xls"),
              row.names = FALSE,quote=FALSE,sep="\t")
  
  write.table(cJun.231.1833.ChipSeq.RNASeq,file=file.path(output.file.dir,"cJun_231_1833_DBA_DGE.xls"),
              row.names = FALSE,quote=FALSE,sep="\t")
  
  write.table(p27.231.1833.ChipSeq.RNASeq,file=file.path(output.file.dir,"p27_231_1833_DBA_DGE.xls"),
              row.names = FALSE,quote=FALSE,sep="\t")
}

addTestFunction4HeatmapChipSeq <- function(GenomicRanges, rtracklayer, IRanges, dataDirectory, chipseq, VennDiagram) {
  library(GenomicRanges)
  library(rtracklayer)
  library(IRanges)
  
  input = import.bed(file.path(dataDirectory, 'ES_input_filtered_ucsc_chr6.bed'),
                     asRangedData=FALSE)
  
  rep1 = import.bed("/Volumes/Bioinformatics$/2017/DannyNewData/PeakCall/2017-03-02-01_S5_R1-MDA-MB-231-DD-1-cJun_hs_1.00e-05_macs142_peaks.bed",genome = "hg19")
                    
  rep2 = import.bed("/Volumes/Bioinformatics$/2017/DannyNewData/PeakCall/2017-03-02-02_S2_R1-MDA-MB-231-DD-2-cJun_hs_1.00e-05_macs142_peaks.bed",genome = "hg19")
  
  library(chipseq)
  estimate.mean.fraglen(rep1)
  
  prepareChIPseq = function(reads){
  
  frag.len = median( estimate.mean.fraglen(reads) )
  cat( paste0( 'Median fragment size for this library is ', round(frag.len)))
  reads.extended = resize(reads, width = frag.len)
  return( trim(reads.extended) )
  }
  
  rep1 = prepareChIPseq(rep1)
  
  ovlp = findOverlaps(rep1,rep2)
  ovlp = as.data.frame(ovlp)
  
  ov = min( length(unique(ovlp$queryHits)), length(unique(ovlp$subjectHits)) )
  
  
  library(VennDiagram)
  draw.pairwise.venn(
    area1=length(rep1),
    area2=length(rep2),
    cross.area=ov,
    category=c("rep1", "rep2"),
    fill=c("steelblue", "blue3"),
    cat.cex=0.7)
}


#input.region.bed.dir = "~/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/ChipSeq/heatmap"
#plotHeatMapUsedeepTools("~/BamCompare","/projects/ctsi/bbc/aimin","hg19_gene.bed",)
#
plotHeatMapUsedeepTools <- function(input.sample.file,input.bw.file.dir,input.region.bed.dir,select.region.bed,output.file.dir){
  
  bw.file.sample.label <- mapBw3Sample(input.sample.file,input.bw.file.dir)
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
#/projects/ctsi/bbc/aimin/annotation/hg19_gene.bed
  
   bed.file.list <- list.files(
        input.region.bed.dir,
        pattern = "*.bed$",
        all.files = TRUE,
        full.names = TRUE,
        recursive = TRUE,
        include.dirs = TRUE)
  
  bw.file.list <- as.character(bw.file.sample.label$file.bw)
  samplesLabel <- as.character(bw.file.sample.label$sampleLabel)
  
  print(bed.file.list)
  
if(select.region.bed!=""){
input.beds = bed.file.list[-grep(paste(select.region.bed,collapse = "|"), bed.file.list)] 
input.beds = paste(input.beds,collapse = " ")
}else{
input.beds = paste(bed.file.list,collapse = " ")
}

  print(input.beds)
  
cmd = "computeMatrix reference-point --referencePoint TSS -b 4000 -a 4000 -R"
input.beds=input.beds
input.bw.file=paste(bw.file.list,collapse = " ")
input.samplesLabel=paste(samplesLabel,collapse = " ")
outFileNameMatrix=file.path(output.file.dir,"matrix1_cJun_p27_TSS.gz")
outFileSortedRegions=file.path(output.file.dir,"regions_cJun_p27_genes.bed")
outHeatMapFile=file.path(output.file.dir,"heatmap_cJun_p27.png")

cmd1=paste(cmd,input.beds,"-S",input.bw.file,"--skipZero","-o",outFileNameMatrix,"--outFileSortedRegions",outFileSortedRegions,"--samplesLabel",input.samplesLabel,collapse = " ")

cmd2=paste("plotHeatmap -m",outFileNameMatrix,"-out",outHeatMapFile,collapse = " ")

# 
# testFiles/genes.bed -S testFiles/log2ratio_H3K4Me3_chr19.bw --skipZeros 
# -o matrix1_H3K4me3_l2r_TSS.gz
# --outFileSortedRegions regions1_H3K4me3_l2r_genes.bed
# 
# deepTools2.0/bin/computeMatrix scale-regions \
#   -R genes_chr19_firstHalf.bed genes_chr19_secondHalf.bed \ # separate multiple files with spaces
#   -S testFiles/log2ratio_*.bw  \ or use the wild card approach
#   -b 3000 -a 3000 \
#   --regionBodyLength 5000 \
#   --skipZeros -o matrix2_multipleBW_l2r_twoGroups_scaled.gz \
#   --outFileNameMatrix matrix2_multipleBW_l2r_twoGroups_scaled.tab \
#   --outFileSortedRegions regions2_multipleBW_l2r_twoGroups_genes.bed
# 
# computeMatrix reference-point \ # choose the mode
#        --referencePoint TSS \ # alternatives: TES, center
#        -b 3000 -a 10000 \ # define the region you are interested in
#        -R testFiles/genes.bed \
#        -S testFiles/log2ratio_H3K4Me3_chr19.bw  \
#        --skipZeros \
#        -o matrix1_H3K4me3_l2r_TSS.gz \ # to be used with plotHeatmap and plotProfile
#        --outFileSortedRegions regions1_H3K4me3_l2r_genes.bed
#cmd.java.2="export LD_LIBRARY_PATH=/usr/lib/jvm/java-1.7.0-openjdk-1.7.0.45.x86_64/jre/lib/amd64/server:$LD_LIBRARY_PATH"
#cmd1=paste(cmd.java.2,cmd1,sep=";")
system(cmd1)

#print(cmd2)

}


#R -e 'library(PathwaySplice);library(ChipSeq);ChipSeq:::submitJob4plotHeatMapUsedeepTools("~/BamCompare","/projects/ctsi/bbc/aimin/annotation","hg19_gene.bed","/scratch/projects/bbc/aiminy_project/Danny_ChipSeq_heatmap")'

submitJob4plotHeatMapUsedeepTools <- function(input.bw.file.dir,input.region.bed.dir,select.region.bed,output.file.dir){
  
  #Sys.setenv(JAVA_HOME='/usr/lib/jvm/java-1.7.0-openjdk-1.7.0.45.x86_64/jre/lib/amd64/server')
  
  if (!dir.exists(output.file.dir))
  {
    dir.create(output.file.dir, recursive = TRUE)
  }
  
  job.name <- "HeapMapPlot"
  
  Rfun1 <- 'library(ChipSeq);re <- ChipSeq:::plotHeatMapUsedeepTools('
  
  Rinput <- paste0('\\"',input.bw.file.dir,'\\",',
                   '\\"',input.region.bed.dir,'\\",',
                   '\\"',select.region.bed,'\\",',
                   '\\"',output.file.dir,'\\"')
  Rfun2 <- ')'
  
  Rfun <-paste0(Rfun1,Rinput,Rfun2)
  
  #cmd.java.1="module load java/1.8.0_60"
  #cmd.java.1="R CMD javareconf -e"
 # cmd.java.2="export LD_LIBRARY_PATH=/usr/lib/jvm/java-1.7.0-openjdk-1.7.0.45.x86_64/jre/lib/amd64/server:$LD_LIBRARY_PATH"
  
  #cmd.java ='export JAVA_HOME="/usr/lib/jvm/java-1.7.0-openjdk-1.7.0.45.x86_64/jre"'
  
  cmd.gff <- PathwaySplice:::createBsubJobArrayRfun(Rfun,job.name,wait.job.name=NULL)
  
  #cmd2=paste(cmd.java.2,cmd.gff,sep=";")
  
  system(cmd.gff)

}

# Welcome to Rattler, *please* read these important system notes:
#   
#   --> Rattler is currently running the SLURM resource manager to
# schedule all compute resources. Example SLURM job scripts are
# available on the system at /cm/shared/docs/slurm/current/job_scripts
# 
# To run an interactive shell, issue:
#   srun -p hipri -t 0:30:00 -n 36 -N 1 --pty /bin/bash -l
# 
# To submit a batch job, issue:       sbatch job_script
# To show all queued jobs, issue:     squeue
# To kill a queued job, issue:        scancel <jobId>
#   
#   See "man slurm" for more detailed information.
# 
# --> Rattler has 3 main queues with the following timelimits:
#   * hipri    (2 days)
# * iq       (1 day)
# * requestq (infinite) - Requires a valid slurm reservation to use
# 
# --> To see which software packages are available issue: "module avail"
# 
# Use the following commands to adjust your environment:
#   'module add <module>'     - adds a module to your environment for this session
# 'module initadd <module>' - configure module to be loaded at every login
# 

#plotHeatMapUsedeepTools("/scratch/projects/bbc/aiminy_project/Danny_ChipSeq_heatmap/matrix1_cJun_p27_TSS.gz","/scratch/projects/bbc/aiminy_project/Danny_ChipSeq_heatmap")
#
usePlotHeatmap <- function(input.matrix.file,output.file.dir){
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  #/projects/ctsi/bbc/aimin/annotation/hg19_gene.bed
  
  outHeatMapFile=file.path(output.file.dir,"heatmap_cJun_p27.png")
  
  cmd2=paste("plotHeatmap -m",input.matrix.file,"-out",outHeatMapFile,collapse = " ")
  
  system(cmd2)
  
}

#R -e 'library(PathwaySplice);library(ChipSeq);ChipSeq:::submitJob4usePlotHeatmap("/scratch/projects/bbc/aiminy_project/Danny_ChipSeq_heatmap/matrix1_cJun_p27_TSS.gz","/scratch/projects/bbc/aiminy_project/Danny_ChipSeq_heatmap")'

submitJob4usePlotHeatmap <- function(input.matrix.file,output.file.dir){
  
  if (!dir.exists(output.file.dir))
  {
    dir.create(output.file.dir, recursive = TRUE)
  }
  
  job.name <- "plotHeatmap"
  
  Rfun1 <- 'library(ChipSeq);re <- ChipSeq:::usePlotHeatmap('
  
  Rinput <- paste0('\\"',input.matrix.file,'\\",',
                   '\\"',output.file.dir,'\\"')
  Rfun2 <- ')'
  
  Rfun <-paste0(Rfun1,Rinput,Rfun2)
  
  cmd.gff <- PathwaySplice:::createBsubJobArrayRfun(Rfun,job.name,wait.job.name=NULL)
  
  system(cmd.gff)
}

# input.bam.file.dir="/Volumes/Bioinformatics$/2017/DannyNewData/BindDiff/common_peaks_bed/Annotation"
# out.dir.name = "~/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/ChipSeq/heatmap"
# generateBed4HeatMap(input.bam.file.dir,out.dir.name)
#
generateBed4HeatMap <- function(input.bam.file.dir,out.dir.name) {
  
  if (!dir.exists(out.dir.name))
  {
    dir.create(out.dir.name, recursive = TRUE)
  }
  
  file.list <- list.files(
    input.bam.file.dir,
    pattern = "*xls$",
    all.files = TRUE,
    full.names = TRUE,
    recursive = TRUE,
    include.dirs = TRUE)
  
  file.table  <- lapply(file.list,function(u){
    y <- basename(u)
    pos <- regexpr("\\_", y)
    pos <- pos - 1
    y <- substr(y, 1, pos)
    x <-fread(u)
    x<- cbind(x,rep(y,dim(x)[1]))
    colnames(x)[dim(x)[2]]="sample_name"
    x
    })
  
  # library(dplyr)
  # x <- data_frame(i = c("a","b","c"), j = 1:3)
  # y <- data_frame(i = c("b","c","d"), k = 4:6)
  # z <- data_frame(i = c("c","d","a"), l = 7:9)
  
  common.gene.peak <- file.table %>%
    Reduce(function(dtf1,dtf2) inner_join(dtf1,dtf2,by="SYMBOL"),.)
  
  x <- data.table::subset(common.gene.peak,select=c("geneChr.x","geneStart.x","geneEnd.x","geneStrand.x","SYMBOL"))
  
  data.table::setkey(x,NULL)
  xx <- unique(x)
  xxx<-xx[complete.cases(xx),]
  
  # gene.231dd.231.1833.1833shp27 <- read.xlsx(
  # "/Volumes/Bioinformatics$/2017/DannyNewData/AnnotationNew4DBA/table_S3_headerless.xlsx",
  # sheetIndex = 1,header = FALSE)
  # colnames(gene.231dd.231.1833.1833shp27)=c("SYMBOL","FC_231DD_231","p_231DD_231","FC_1833_231","p_1833_231","FC_1833shp27_1833","p_1833shp27_1833")
  #save(gene.231dd.231.1833.1833shp27,file="./data/gene.RData")
  gene.45 <- merge(gene.231dd.231.1833.1833shp27,xxx,by="SYMBOL",sort = FALSE)
  
  gene.45.1 <- gene.45[,c(8,9,10,11,1)]
  colnames(gene.45.1)=c("chr","start","end","strand","SYMBOL")
  
  gene.45.1[which(gene.45.1$strand==1),]$strand="+"
  gene.45.1[which(gene.45.1$strand==2),]$strand="-"
  
  
  gene.45.1.sorted.by.chr <- gene.45.1[order(gene.45.1$chr),]
  
  gene.45.1$chr<-paste0("chr",gene.45.1$chr)
  gene.45.1.sorted.by.chr$chr<-paste0("chr",gene.45.1.sorted.by.chr$chr)
  
  write.table(gene.45.1,file=file.path(out.dir.name,"gene_200.bed"),row.names = FALSE,col.names = FALSE,quote=FALSE,sep="\t")
  write.table(gene.45.1.sorted.by.chr,file=file.path(out.dir.name,"gene_200_sorted.bed"),row.names = FALSE,col.names = FALSE,quote=FALSE,sep="\t")
}


#R -e 'library(PathwaySplice);library(ChipSeq);ChipSeq:::bashJob4plotHeatMapUsedeepTools("~/SampleID_INFO_ChIP_new_Danny.csv","~/BamCompare","~/ChipSeqBed",select.region.bed=NULL,"/scratch/projects/bbc/aiminy_project/Danny_ChipSeq_heatmapAddlabel")'

#R -e 'library(PathwaySplice);library(ChipSeq);ChipSeq:::bashJob4plotHeatMapUsedeepTools("~/SampleID_INFO_ChIP_new_Danny.csv","~/BamCoverage","~/ChipSeqBed_p27_cJun",select.region.bed=NULL,"/scratch/projects/bbc/aiminy_project/Danny_ChipSeq_heatmap4p27_cJun")'


#R -e 'library(PathwaySplice);library(ChipSeq);ChipSeq:::bashJob4plotHeatMapUsedeepTools("/nethome/axy148/R/lib64/R/library/ChipSeq/extdata/zhao_data.csv","~/BamCompareZhao","~/cJun_gene_bed",select.region.bed=NULL,"/scratch/projects/bbc/aiminy_project/zhao_ChipSeq_heatmap4cJun")'

bashJob4plotHeatMapUsedeepTools <- function(input.sample.file,input.bw.file.dir,input.region.bed.dir,select.region.bed,output.file.dir){
  
  #Sys.setenv(JAVA_HOME='/usr/lib/jvm/java-1.7.0-openjdk-1.7.0.45.x86_64/jre/lib/amd64/server')
  
  if (!dir.exists(output.file.dir))
  {
    dir.create(output.file.dir, recursive = TRUE)
  }
  
  job.name <- "generateMatrix"
  
  Rfun1 <- 'library(ChipSeq);re <- ChipSeq:::plotHeatMapUsedeepTools('
  
  
  Rinput <- paste0('\\"',input.sample.file,'\\",',
                     '\\"',input.bw.file.dir,'\\",',
                    '\\"',input.region.bed.dir,'\\",',
                   '\\"',select.region.bed,'\\",',
                   '\\"',output.file.dir,'\\"')

  Rfun2 <- ')'
  
  Rfun <-paste0(Rfun1,Rinput,Rfun2)
  
  #cmd.java.1="module load java/1.8.0_60"
  #cmd.java.1="R CMD javareconf -e"
  # cmd.java.2="export LD_LIBRARY_PATH=/usr/lib/jvm/java-1.7.0-openjdk-1.7.0.45.x86_64/jre/lib/amd64/server:$LD_LIBRARY_PATH"
  
  #cmd.java ='export JAVA_HOME="/usr/lib/jvm/java-1.7.0-openjdk-1.7.0.45.x86_64/jre"'
  
  cmd.gff <- PathwaySplice:::createBsubJobArrayRfun(Rfun,job.name,wait.job.name=NULL)
  
  #cmd2=paste(cmd.java.2,cmd.gff,sep=";")
  
  system(cmd.gff)
  
  job.name <- "plotHeatmap"
  
  Rfun1 <- 'library(ChipSeq);re <- ChipSeq:::usePlotHeatmap('
  input.matrix.file <- file.path(output.file.dir,"matrix1_cJun_p27_TSS.gz")
  
  Rinput <- paste0('\\"',input.matrix.file,'\\",',
                   '\\"',output.file.dir,'\\"')
  Rfun2 <- ')'
  
  Rfun <-paste0(Rfun1,Rinput,Rfun2)
  
  cmd.gff <- PathwaySplice:::createBsubJobArrayRfun(Rfun,job.name,wait.job.name="generateMatrix")
  
  system(cmd.gff)
  
}

#input.bw.dir="~/BamCompare"
#input.sample.file="~/SampleID_INFO_ChIP_new_Danny.csv"
#
#mapBw3Sample(input.sample.file,input.bw.dir)
#
mapBw3Sample <- function(input.sample.file,input.bw.dir) {
  
  file.1 <- list.files(input.bw.dir,pattern=".bw$",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
  
  sample.file <- fread(input.sample.file)
  
  file.2<-cbind(unlist(lapply(file.1,function(u){x<-basename(u);p1 <- regexpr("\\_", x);p2 <- regexpr("\\.", x);xx <- substr(x,p1+1,p2-1)})),file.1)
  
  colnames(file.2)=c("ID","file.bw")
  file.3 <- merge(file.2,sample.file,by="ID",sort=F)
  sampleLabel= paste(gsub(" ", "-", file.3$Type_Cell), file.3$Type_TF, sep = "-")
  sampleLabel=gsub("MDA-MB-","",sampleLabel)
  file.4 <- cbind(file.3,sampleLabel)
  file.5 <- file.4[,c(2,6)]
  file.5
}

makeHeatMapByme <- function() {
  compute.matrix.file="/Users/axy148/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/ChipSeq/heatmap/Danny_ChipSeq_heatmapBasedSortedBed/matrix1_cJun_p27_TSS.gz"
  
  table = read.table(gzfile(compute.matrix.file), skip=1)
  matrix = table[7:dim(table)[2]]
  image(1:dim(matrix)[2], 1:dim(matrix)[1], t(as.matrix(matrix)),  axes=FALSE, xlab='sample', ylab='gene')
  axis(2, at=1:length(table$V1), labels=table$V4, las = 1)
}

# input.bam.file.dir="/Volumes/Bioinformatics$/2017/DannyNewData/BindDiff/common_peaks_bed/Annotation"
# out.dir.name = "~/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/ChipSeq/heatmap"
# gene.full <- generateBed4HeatMap2(input.bam.file.dir,out.dir.name)
#
# input.bam.file.dir="~/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/ChipSeq/heatmap/Peak2Bed/p27"
# out.dir.name = "~/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/ChipSeq/heatmap/Peak2Bed"
# gene.full <- generateBed4HeatMap2(input.bam.file.dir,out.dir.name,"p27")
#
# input.bam.file.dir="~/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/ChipSeq/heatmap/Peak2Bed/cJun"
# out.dir.name = "~/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/ChipSeq/heatmap/Peak2Bed"
# gene.full <- generateBed4HeatMap2(input.bam.file.dir,out.dir.name,"cJun")
#

generateBed4HeatMap2 <- function(input.bam.file.dir,out.dir.name,Ab) {
  
  out.dir.name <-file.path(out.dir.name,paste0(Ab,"_gene_bed"))
  
  if (!dir.exists(out.dir.name))
  {
    dir.create(out.dir.name, recursive = TRUE)
  }
  
  file.list <- list.files(
    input.bam.file.dir,
    pattern = "*xls$",
    all.files = TRUE,
    full.names = TRUE,
    recursive = TRUE,
    include.dirs = TRUE)
  
  file.table  <- lapply(file.list,function(u){
    y <- basename(u)
    pos <- regexpr("\\_", y)
    pos <- pos - 1
    y <- substr(y, 1, pos)
    x <-fread(u)
    x<- cbind(x,rep(y,dim(x)[1]))
    colnames(x)[dim(x)[2]]="sample_name"
    x
  })
  
  # library(dplyr)
  # x <- data_frame(i = c("a","b","c"), j = 1:3)
  # y <- data_frame(i = c("b","c","d"), k = 4:6)
  # z <- data_frame(i = c("c","d","a"), l = 7:9)
  
  common.gene.peak <- file.table %>%
    Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="SYMBOL"),.)
  
  x <- subset(common.gene.peak,select=c("geneChr.x","geneStart.x","geneEnd.x","geneStrand.x","SYMBOL"))
  
  setkey(x,NULL)
  xx <- unique(x)
  xxx<-xx[complete.cases(xx),]
  
  # gene.231dd.231.1833.1833shp27 <- read.xlsx(
  # "/Volumes/Bioinformatics$/2017/DannyNewData/AnnotationNew4DBA/table_S3_headerless.xlsx",
  # sheetIndex = 1,header = FALSE)
  # colnames(gene.231dd.231.1833.1833shp27)=c("SYMBOL","FC_231DD_231","p_231DD_231","FC_1833_231","p_1833_231","FC_1833shp27_1833","p_1833shp27_1833")
  #save(gene.231dd.231.1833.1833shp27,file="./data/gene.RData")
  #gene.45 <- merge(gene.231dd.231.1833.1833shp27,xxx,by="SYMBOL",sort = FALSE)
  
  gene.45 <- xxx
  gene.45.1 <- gene.45
  colnames(gene.45.1)=c("chr","start","end","strand","SYMBOL")

  gene.45.1[which(gene.45.1$strand==1),]$strand="+"
  gene.45.1[which(gene.45.1$strand==2),]$strand="-"


  gene.45.1.sorted.by.chr <- gene.45.1[order(gene.45.1$chr),]

  gene.45.1$chr<-paste0("chr",gene.45.1$chr)
  gene.45.1.sorted.by.chr$chr<-paste0("chr",gene.45.1.sorted.by.chr$chr)

  write.table(gene.45.1,file=file.path(out.dir.name,paste0("gene_",Ab,".bed")),row.names = FALSE,col.names = FALSE,quote=FALSE,sep="\t")
  write.table(gene.45.1.sorted.by.chr,file=file.path(out.dir.name,paste0("gene_",Ab,"_sorted.bed")),row.names = FALSE,col.names = FALSE,quote=FALSE,sep="\t")
  
  # > length(unique(file.table[[1]]$SYMBOL))
  # [1] 417
  # > length(unique(file.table[[2]]$SYMBOL))
  # [1] 222
  # > length(unique(file.table[[3]]$SYMBOL))
  # [1] 458
  # 
  # 
  # unique(c(unique(file.table[[1]]$SYMBOL),unique(file.table[[2]]$SYMBOL),unique(file.table[[3]]$SYMBOL)))
  gene.45.1.sorted.by.chr
}



# /deepTools-1.5/bin/bamCompare --bamfile1 ChIP.bam --bamfile2 Input.bam \
# --binSize 25 --fragmentLength 200 --missingDataAsZero no \
# --ratio log2 --scaleFactorsMethod SES -o log2ratio_ChIP_vs_Input.bw

#' bsub -P bbc -J "BamCoverage" -o %J.BamCoverage.log -e %J.BamCoverage.err -W 72:00 -n 32 -q parallel -R 'rusage[mem= 16000 ] span[ptile= 16 ]' -u aimin.yan@med.miami.edu R -e 'library(ChipSeq);re <- ChipSeq:::useBamCoverage("/scratch/projects/bbc/aiminy_project/DannyNewData2/SampleID_INFO_ChIP_new_Danny.csv","/nethome/axy148/_visualization/sorted_bam_files_2.txt","~/BamCoverage")'

useBamCoverage <- function(input.sample.file,input.bam.file,output.dir)
{
  
  re <- GetSampleInfo(input.sample.file, input.bam.file)
  
  cellInfo <- re$y
  
  # output.dir.name = dirname(input.sample.file)
  
  if (!dir.exists(output.dir))
  {
    dir.create(output.dir, recursive = TRUE)
  }
  
  temp3 = output.dir
  
  # cmd9 = 'ngs.plot.r -G' cmd10 = '-R' cmd11 = '-C' cmd12 = '-O' cmd13 = '-T'
  # cmd14 = '-L' cmd15 = '-RR' cmd16 = '-CD' cmd17= '-GO'
  
  
  
  cellInfo.run <- lapply(1:length(cellInfo), function(u,cellInfo, 
                                                      temp3)
  {
    
    x.name = cellInfo[[u]]$name
    
    es <- cellInfo[[u]]$es
    
    x.input <- es[es$Type_TF == "Input", ]$file.name
    
    x.sample <- es[es$Type_TF != "Input", ]
    
    #print(x.sample)
    #print(x.input)
    
    x.run <- apply(x.sample, 1, function(x, x.input, temp3)
    {
      
      y <- x
      
      ID <- y[1]
      Type_Cell <- y[2]
      Type_TF <- y[3]
      Cell_TF <- y[4]
      file.name <- y[5]
      xx <- file.name
      xx.name = paste(ID, gsub(" ", "-", Type_Cell), Type_TF, sep = "-")
      
      # ~/python/Python-2.7.11/python  ~/NGS_tools/deepTools/bin/bamCompare  --bamfile1 /scratch/projects/bbc/aiminy_project/DannyNewNgsPlot/2017-03-02-03_S11_R1.marked_sorted.bam --bamfile2 /scratch/projects/bbc/aiminy_project/DannyNewNgsPlot/2017-03-02-17_S1_R1.marked_sorted.bam --binSize 25 --ratio log2 -o ~/BamCompare/log2ratio_2017-03-02-03_S11_R1.marked_sorted.bam_vs_2017-03-02-17_S1_R1.marked_sorted.bam.bw 
      
      #bamCoverage -b reads.bam -o coverage.bw
      
      #bamCoverage --bam a.bam -o a.SeqDepthNorm.bw \
      #--binSize 10
      #--normalizeTo1x 2150570000
      #--ignoreForNormalization chrX
      #--extendReads
      
      cmd1 <- paste("~/python/Python-2.7.11/python  ~/NGS_tools/deepTools/bin/bamCoverage --bam",xx,sep=" ")
      cmd2 <- "--binSize 25 --normalizeTo1x 2451960000 --ignoreForNormalization chrX"
      cmd3 <- "-o"
      
      cmd4 <- paste(cmd1,cmd2,cmd3,sep=" ") 
      
      cmd5 <- file.path(output.dir,paste0("To1x_",basename(as.character(xx)),"_vs_","no_input",".bw"))
      
      
      cmd6 <- paste(cmd4,cmd5,sep=" ") 
      
      cmd6
      
      #cat(cmd6, "\n")
      #cat("\n")
    }, x.input, temp3)
    
    x.run
    
  }, cellInfo,temp3)
  
  
  names(cellInfo.run) = unlist(lapply(cellInfo, function(u)
  {
    u$name
  }))
  
  zzz <- unlist(cellInfo.run)
  
  lapply(1:length(zzz), function(u, zzz)
  {
    
    cat(as.character(zzz[u][[1]]), "\n")
    cat("\n")
    
    system(as.character(zzz[u][[1]]))
    
  }, zzz)
  # # dir.name=temp3 dir.name=reformatPath(dir.name)
  # 
  # file.name = file.path(temp3, dir(temp3, recursive = TRUE))
  # 
  # file.name.2 <- as.list(file.name)
  # 
  # 
  # # names(file.name.2) = unlist(lapply(file.name.2, function(u) { u$name }))
  # 
  # zzz <- unlist(file.name.2)
  # 
  # lapply(1:length(zzz), function(u, zzz)
  # {
  #   
  #   dir.name = dirname(zzz[u][[1]])
  #   file_name = file_path_sans_ext(basename(zzz[u][[1]]))
  #   
  #   cmd = paste("ngs.plot.r -G hg19 -R tss -C", zzz[u][[1]], "-O", file.path(dir.name, 
  #                                                                            paste0(file_name, ".tss")), "-T", file_name, "-L 4000 -RR 1 -CD 1 -GO total", 
  #               sep = " ")
  #   
  #   
  #   
  #   # system(as.character(zzz[u][[1]])) job.name = paste0('bamPlot.', u)
  #   # cmd.pegasus = usePegasus(job.option, Wall.time = '72:00',cores = 32,Memory
  #   # = 16000,span.ptile = 16,job.name) cmd2 = paste(cmd.pegasus,cmd,sep = ' ')
  #   
  #   cmd2 = cmd
  #   cat(cmd2, "\n")
  #   cat("\n")
  #   system(cmd2)
  # }, zzz)
  # 
  # 
  # re <- list(cellInforun = cellInfo.run, zzz = zzz)
  # 
  # # AnntationUsingChipSeeker(temp3, 'peaks.bed', temp3, DD = 5000)
  
  #return(re)
  
}

#' output.file.dir="~/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/ChipSeq/heatmap/Peak2Bed"
#' #mcf7=mcf2
#' re<-getTargetGene4Ab(mcf2,"yes","p27",output.file.dir)
#' re<-getTargetGene4Ab(mcf2,"yes","cJun",output.file.dir)
#' 
getTargetGene4Ab<-function(mcf7,Mergereplicates=c("yes","no"),Ab,output.file.dir){
  
  output.file.dir = file.path(output.file.dir,Ab)
  if (!dir.exists(output.file.dir))
  {
    dir.create(output.file.dir,recursive = TRUE)
  }
  
  temp<-mcf7
  
  sampID.v<-colnames(temp$class) 
  sampID.v.2<-unlist(lapply(1:length(sampID.v),function(u,sampID.v){
    
    x=sampID.v[u]
    
    y=x
    y
  },sampID.v))
  
  
  colnames(temp$class)<-sampID.v.2
  temp$class[1,]<-sampID.v.2
  
  #temp<-dba(temp,mask =c(1,2,13,14,11,12,15,16))
  
  #Merge the replicates of each set 
  if(Mergereplicates=="yes"){
    temp2<-dba.peakset(temp,consensus = -DBA_REPLICATE)
    #temp22<-dba(temp2,mask =c(9,16,17:23))
    print(temp2)
    ##need to set the flexiable number identified from data sets
    t <- temp2
    temp22<-dba(t,mask =which(t$class[which(rownames(t$class)=="Replicate"),]=="1-2"))
  }else{
    temp22<-temp
  }
  
  print(temp22)
  #temp22 <- re
  
   Tissue1<-temp22$class[2,]
   Tissue2<-unique(temp22$class[2,])
   
   TF<-unique(temp22$class[3,])
   TF.n<-length(TF)
    
   for(i in 1:length(Tissue2)){
     
     #for(j in 1:length(TF)){
     p <- which(temp22$class[2,]==Tissue2[i]&temp22$class[3,]==Ab)
     cat(Tissue2[i]," ",Ab," ",p,"\n")
     if(length(p)!=0){
     #as.data.frame(dba.peakset(temp22,peaks=1,bRetrieve=TRUE))
     dba.peakset(temp22,peaks=as.integer(p),bRetrieve=TRUE, writeFile = file.path(output.file.dir,paste0(Tissue2[i],"_",Ab,".bed")))
    }
     
   }
   
   AnntationUsingChipSeeker3(output.file.dir,".bed",output.file.dir
                            ,txdb="hg19",DD=5000)
   
#   AnntationUsingChipSeeker2
  # temp3=file.path(output.file.dir,"Venn")
  # 
  # if(!dir.exists(temp3))
  # {
  #   dir.create(temp3)
  # }
  # 
  # for(i in 1:length(Tissue2)){
  #   po<-which(Tissue1 %in% Tissue2[i])
  #   
  #   print(po)
  #   
  #   if(length(po)==2)
  #   {
  #     png(file.path(temp3,paste0(paste0(colnames(temp22$class)[po],collapse = "-vs-"),".png")))
  #     dba.plotVenn(temp22,mask=po,main=Tissue2[i])
  #     dev.off()
  #     
  #   }else if(length(po)==4){
  #     po1<-po[c(1,3)]
  #     po2<-po[c(2,4)]
  #     
  #     png(file.path(temp3,paste0(paste0(colnames(temp22$class)[po1],collapse = "-vs-"),".png")))
  #     dba.plotVenn(temp22,mask=po1,main=Tissue2[i])
  #     dev.off()
  #     
  #     png(file.path(temp3,paste0(paste0(colnames(temp22$class)[po2],collapse = "-vs-"),".png")))
  #     dba.plotVenn(temp22,mask=po2,main=Tissue2[i])
  #     dev.off()
  #     
  #   }else
  #   {
  #     cat(paste0("For ",Tissue2[i],": Only one peak profile for ",TF.n," TFs\n"))
  #   }
  # }
  # 
  # p.common<-lapply(1:length(Tissue2),function(u,Tissue1,Tissue2,temp22){
  #   
  #   po<-which(Tissue1 %in% Tissue2[u])
  #   
  #   if(length(po)==2)
  #   {
  #     common.peaks<-dba.overlap(temp22,mask=po)
  #     y<-common.peaks$inAll
  #   }else if(length(po)==4){
  #     po1<-po[c(1,3)]
  #     po2<-po[c(2,4)]
  #     
  #     common.peaks.1<-dba.overlap(temp22,mask=po1)
  #     y1<-common.peaks.1$inAll
  #     
  #     common.peaks.2<-dba.overlap(temp22,mask=po2)
  #     y2<-common.peaks.2$inAll
  #     
  #     y<-list(y1=y1,y2=y2)
  #     
  #     names(y)[1]<-paste0(colnames(temp22$class)[po1],collapse = "-vs-")
  #     names(y)[2]<-paste0(colnames(temp22$class)[po2],collapse = "-vs-")
  #     
  #   }else{
  #     y<-NULL}
  #   y
  # },Tissue1,Tissue2,temp22)
  # 
  # names(p.common)<-Tissue2
  # 
  # p.common<-unlist(p.common,recursive = F)
  # 
  # p.common<-p.common[lapply(p.common,length) > 0]
  # 
  # if(!dir.exists(file.path(output.file.dir,"common_peaks_bed")))
  # {
  #   dir.create(file.path(output.file.dir,"common_peaks_bed"))
  #   dir.create(file.path(output.file.dir,"common_peaks_bed","ucsc"))
  #   dir.create(file.path(output.file.dir,"common_peaks_bed","igv"))
  # }
  # 
  # #output common peaks to bed files
  # 
  # lapply(1:length(p.common),function(u,p.common,output.file.dir){
  #   x=p.common[[u]]
  #   
  #   x_name=names(p.common)[u]
  #   
  #   df <- data.frame(seqnames=seqnames(x),
  #                    #starts=start(x)-1,
  #                    starts=start(x),
  #                    ends=end(x),
  #                    names=c(rep(".", length(x))),
  #                    scores=elementMetadata(x)[,1],
  #                    strands=strand(x))
  #   
  #   #assign strand
  #   df.str <- data.frame(seqnames=seqnames(x),
  #                        #starts=start(x)-1,
  #                        starts=start(x),
  #                        ends=end(x),
  #                        names=c(rep(".", length(x))),
  #                        scores=elementMetadata(x)[,1],
  #                        strands=c(rep(".", length(x))))
  #   
  #   df.str.1<-df.str[-grep("random",df.str$seqnames),]
  #   
  #   df.str.2<-df.str.1
  #   
  #   df.str.3<-df.str.2[-grep("chrUn",df.str.2$seqnames),]
  #   
  #   write.table(df,file=file.path(output.file.dir,"common_peaks_bed",paste0(x_name,"_cp_with_header.bed")),
  #               col.names=TRUE,row.names = FALSE,quote=FALSE,sep="\t")
  #   
  #   write.table(df,file=file.path(output.file.dir,"common_peaks_bed","ucsc",paste0(x_name,"_4_ucsc.bed")),
  #               col.names=FALSE,row.names = FALSE,quote=FALSE,sep="\t")
  #   
  #   write.table(df,file=file.path(output.file.dir,"common_peaks_bed",paste0(x_name,"_common_peaks.bed")),
  #               col.names=FALSE,row.names = FALSE,quote=FALSE,sep="\t")
  #   
  #   write.table(df,file=file.path(output.file.dir,"common_peaks_bed","igv",paste0(x_name,"_4_igv.bed")),
  #               col.names=FALSE,row.names = FALSE,quote=FALSE,sep="\t")
  #   
  # },p.common,output.file.dir)
  # 
  # AnntationUsingChipSeeker(file.path(output.file.dir,"common_peaks_bed","igv"),"bed",file.path(output.file.dir,"common_peaks_bed")
  #                          ,txdb="hg19",DD=5000,distanceToTSS_cutoff=10000)
  # 
  # return(p.common)
}

#bsub -P bbc -J "RunSppR" -o %J.RunSppR.log -e %J.RunSppR.err -W 72:00 -n 32 -q parallel -R 'rusage[mem= 16000 ] span[ptile= 16 ]' -u aimin.yan@med.miami.edu R -e 'library(ChipSeq);re <- ChipSeq:::useRunSppR("~/SampleID_INFO_ChIP_new_Danny.csv","~/sort_bam","/scratch/projects/bbc/aiminy_project/Danny_ChipSeq_data_QC")'
#
#R/lib64/R/library/ChipSeq/extdata/zhao_data.csv 
#
#bsub -P bbc -J "RunSppR" -o %J.RunSppR.log -e %J.RunSppR.err -W 72:00 -n 32 -q parallel -R 'rusage[mem= 16000 ] span[ptile= 16 ]' -u aimin.yan@med.miami.edu R -e 'library(ChipSeq);re <- ChipSeq:::useRunSppR("R/lib64/R/library/ChipSeq/extdata/zhao_data.csv","~/sort_bam","/scratch/projects/bbc/aiminy_project/zhao_ChipSeq_data_QC")'
#
#bsub -P bbc -J "RunSppR" -o %J.RunSppR.log -e %J.RunSppR.err -W 72:00 -n 32 -q parallel -R 'rusage[mem= 16000 ] span[ptile= 16 ]' -u aimin.yan@med.miami.edu R -e 'library(ChipSeq);re <- ChipSeq:::useRunSppR("~/SampleID_INFO_ChIP_new_Danny.csv","/projects/scratch/bbc/Project/Danny_chip2/BigWig","/scratch/projects/bbc/aiminy_project/Danny_ChipSeq_data_QC_2")'
#
#
#
useRunSppR <- function(input.sample.file,input.bam.file.dir,output.file.dir){
  
  #bam.file.sample.label <- mapBam2Sample2(input.sample.file,input.bam.file.dir)
  
  bam.file.sample.label <- mapBam2Sample3(input.sample.file,input.bam.file.dir)
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  #/projects/ctsi/bbc/aimin/annotation/hg19_gene.bed
  
  bam.file.list <- as.character(bam.file.sample.label$file.bam)
  samplesLabel <- as.character(bam.file.sample.label$sampleLabel)
  
  lapply(1:length(bam.file.list),function(u,bam.file.list,samplesLabel,output.file.dir){
    
    cmd = paste("Rscript ~/phantompeakqualtools/run_spp.R",
                paste0("-c=",bam.file.list[u]),
                paste0("-savp=",file.path(output.file.dir,paste0(samplesLabel[u],"_QC.pdf"))),
                paste0("-out=",file.path(output.file.dir,paste0(samplesLabel[u],"_QC_metrix.txt"))),sep=" ")
    
    system(cmd)
    
  },bam.file.list,samplesLabel,output.file.dir)

}

#R -e 'library(PathwaySplice);library(ChipSeq);ChipSeq:::submitJob4plotHeatMapUsedeepTools("~/BamCompare","/projects/ctsi/bbc/aimin/annotation","hg19_gene.bed","/scratch/projects/bbc/aiminy_project/Danny_ChipSeq_heatmap")'

submitJob4useRunSppR <- function(input.sample.file,input.bam.file.dir,output.file.dir){
  
  #Sys.setenv(JAVA_HOME='/usr/lib/jvm/java-1.7.0-openjdk-1.7.0.45.x86_64/jre/lib/amd64/server')
  
  if (!dir.exists(output.file.dir))
  {
    dir.create(output.file.dir, recursive = TRUE)
  }
  
  job.name <- "QC"
  
  Rfun1 <- 'library(ChipSeq);re <- ChipSeq:::useRunSppR('
  
  Rinput <- paste0('\\"',input.sample.file,'\\",',
                   '\\"',input.bam.file.dir,'\\",',
                   '\\"',output.file.dir,'\\"')
  Rfun2 <- ')'
  
  Rfun <-paste0(Rfun1,Rinput,Rfun2)
  
  cmd.gff <- DoGs:::createBsubJobArrayRfun(Rfun,job.name,wait.job.name=NULL)
  
  system(cmd.gff)
  
}

#input.bw.dir="~/BamCompare"
#input.sample.file="~/SampleID_INFO_ChIP_new_Danny.csv"
#
#mapBw3Sample(input.sample.file,input.bw.dir)
#
mapBam2Sample <- function(input.sample.file,input.bam.dir) {
  
  file.1 <- list.files(input.bam.dir,pattern=".bam$",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
  
  sample.file <- fread(input.sample.file)
  
  file.2<-cbind(unlist(lapply(file.1,function(u){x<-basename(u);p1 <- regexpr("\\_", x);p2 <- regexpr("\\.", x);xx <- substr(x,p1+1,p2-1)})),file.1)
  
  colnames(file.2)=c("ID","file.bam")
  file.3 <- merge(file.2,sample.file,by="ID",sort=F)
  sampleLabel= paste(gsub(" ", "-", file.3$Type_Cell), file.3$Type_TF, sep = "-")
  sampleLabel=gsub("MDA-MB-","",sampleLabel)
  file.4 <- cbind(file.3,sampleLabel)
  file.5 <- file.4[,c(2,6)]
  file.5
}

mapBam2Sample2 <- function(input.sample.file,input.bam.dir) {
  
  file.1 <- list.files(input.bam.dir,pattern=".bam$",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
  
  sample.file <- fread(input.sample.file)
  
  file.2<-cbind(unlist(lapply(file.1,function(u){x<-basename(u);p2 <- regexpr("\\.", x);xx <- substr(x,1,p2-1)})),file.1)
  
  colnames(file.2)=c("ID","file.bam")
  file.3 <- merge(file.2,sample.file,by="ID",sort=F)
  sampleLabel= paste(gsub(" ", "-", file.3$Type_Cell), file.3$Type_TF, sep = "-")
  sampleLabel=gsub("MDA-MB-","",sampleLabel)
  file.4 <- cbind(file.3,sampleLabel)
  file.5 <- file.4[,c(2,6)]
  file.5
}

mapBam2Sample3 <- function(input.sample.file,input.bam.dir) {
  
  file.1 <- list.files(input.bam.dir,pattern=".bam$",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
  
  sample.file <- fread(input.sample.file)
  
  file.2<-cbind(unlist(lapply(file.1,function(u){x<-basename(u);xx <- tools::file_path_sans_ext(x)})),file.1)
  
  colnames(file.2)=c("ID","file.bam")
  file.3 <- merge(file.2,sample.file,by="ID",sort=F)
  sampleLabel= paste(gsub(" ", "-", file.3$Type_Cell), file.3$Type_TF, sep = "-")
  sampleLabel=gsub("MDA-MB-","",sampleLabel)
  file.4 <- cbind(file.3,sampleLabel)
  file.5 <- file.4[,c(2,6)]
  file.5
}

#' dir.name="~/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/ChipSeq/Bba2Bed"
#' input.file.pattern=".bed"
#' out.dir.name="~/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/ChipSeq/AnnotationNew4DBA"
#' txdb="hg19"
#' DD=5000
#' 
#' AnntationUsingChipSeeker3(dir.name,input.file.pattern,out.dir.name,txdb=txdb,DD,distanceToTSS_cutoff=5000, AP=c("Promoter","Intron"))
#'
#' res.promoter <- AnntationUsingChipSeeker3(dir.name,input.file.pattern,out.dir.name,txdb=txdb,DD,distanceToTSS_cutoff=5000,AP=c("Promoter"))
#' 
#' AnntationUsingChipSeeker3(dir.name,input.file.pattern,out.dir.name,txdb=txdb,DD,distanceToTSS_cutoff=5000,AP=c("Intron"))
#' 

AnotationUsingChipSeeker3 <- function(dir.name,input.file.pattern,out.dir.name,txdb=c("hg19","hg38"),DD,distanceToTSS_cutoff=NULL,assignGenomicAnnotation=TRUE,AP=c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic")) {
  
  re<-ParserReadFiles(dir.name,input.file.pattern)
  
  re.bed<-re$input

  re.peaks.only.bed.2 <- re.bed

  txdb<-match.arg(txdb)
  
  switch (txdb,
          hg38 = {
            cat("Use hg38\n")
            txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
          },
          {
            cat("Use hg19\n") 
            txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
          }
  )
  
  APpath <- paste(AP,collapse = "_")
  
  temp3=file.path(out.dir.name,APpath)
  
  if(!dir.exists(temp3)){dir.create(temp3,recursive = TRUE)}
  
  d=DD
  
  peaks.anno.list <- lapply(1:length(re.peaks.only.bed.2),function(u,re.peaks.only.bed.2,d){
    
    peaks=readPeakFile(re.peaks.only.bed.2[[u]],as="data.frame")
    
    print(head(peaks))
    
    if(dim(peaks)[1]>0){
    
    peakAnno <- annotatePeak(re.peaks.only.bed.2[[u]], tssRegion=c(-d, d),
                             TxDb=txdb,assignGenomicAnnotation=assignGenomicAnnotation,genomicAnnotationPriority=AP,annoDb="org.Hs.eg.db")
    
    dropAnnoM <- function (csAnno, distanceToTSS_cutoff) 
    {
      idx <- which(abs(mcols(csAnno@anno)[["distanceToTSS"]]) < 
                     distanceToTSS_cutoff)
      csAnno@anno <- csAnno@anno[idx]
      csAnno@peakNum <- length(idx)
      if (csAnno@hasGenomicAnnotation) {
        csAnno@annoStat <- ChIPseeker:::getGenomicAnnoStat(csAnno@anno)
        csAnno@detailGenomicAnnotation = csAnno@detailGenomicAnnotation[idx, 
                                                                        ]
      }
      csAnno
    }
    
    if(!is.null(distanceToTSS_cutoff)){
    peakAnno <- dropAnnoM(peakAnno,distanceToTSS_cutoff = distanceToTSS_cutoff)
    }else
    {
    peakAnno <- peakAnno 
    }
    
    x_name=names(re.peaks.only.bed.2)[u]
    x_name=tools::file_path_sans_ext(x_name)
    x_name=tools::file_path_sans_ext(x_name)
    
    cat(x_name,"\n")
    png(file.path(temp3,paste0(x_name,"_",d,".png")))
    plotAnnoPie(peakAnno)
    dev.off()
    
    peaks.anno=as.data.frame(peakAnno)
    
    print(head(peaks.anno))
    
    print(colnames(peaks.anno))
    write.table(peaks.anno,file=file.path(temp3,paste0(x_name,"_",d,".xls")),
                row.names = FALSE,quote=FALSE,sep="\t")

    peaks.anno
    }
    
  },re.peaks.only.bed.2,d)
  
}

# R -e 'library(PathwaySplice);library(ChipSeq);ChipSeq:::bashJob4generateBed4HeatMap2("~/cJun","~/","cJun")'

bashJob4generateBed4HeatMap2 <- function(input.bam.file.dir,out.dir.name,Ab){
  
  #Sys.setenv(JAVA_HOME='/usr/lib/jvm/java-1.7.0-openjdk-1.7.0.45.x86_64/jre/lib/amd64/server')
  
  if (!dir.exists(out.dir.name))
  {
    dir.create(out.dir.name, recursive = TRUE)
  }
  
  job.name <- "generateGeneBed"
  
  Rfun1 <- 'library(ChipSeq);re <- ChipSeq:::generateBed4HeatMap2('
  
  
  Rinput <- paste0('\\"',input.bam.file.dir,'\\",',
                   '\\"',out.dir.name,'\\",',
                   '\\"',Ab,'\\"')

  Rfun2 <- ')'
  
  Rfun <-paste0(Rfun1,Rinput,Rfun2)
  
  #cmd.java.1="module load java/1.8.0_60"
  #cmd.java.1="R CMD javareconf -e"
  # cmd.java.2="export LD_LIBRARY_PATH=/usr/lib/jvm/java-1.7.0-openjdk-1.7.0.45.x86_64/jre/lib/amd64/server:$LD_LIBRARY_PATH"
  
  #cmd.java ='export JAVA_HOME="/usr/lib/jvm/java-1.7.0-openjdk-1.7.0.45.x86_64/jre"'
  
  cmd.gff <- PathwaySplice:::createBsubJobArrayRfun(Rfun,job.name,wait.job.name=NULL)
  
  #cmd2=paste(cmd.java.2,cmd.gff,sep=";")
  
  system(cmd.gff)
  
}

# /deepTools-1.5/bin/bamCompare --bamfile1 ChIP.bam --bamfile2 Input.bam \
# --binSize 25 --fragmentLength 200 --missingDataAsZero no \
# --ratio log2 --scaleFactorsMethod SES -o log2ratio_ChIP_vs_Input.bw

#' bsub -P bbc -J "bamCompare" -o %J.bamCompare.log -e %J.Compare.err -W 72:00 -n 32 -q parallel -R 'rusage[mem= 16000 ] span[ptile= 16 ]' -u aimin.yan@med.miami.edu R -e 'library(ChipSeq);re <- ChipSeq:::useBamCompare2("/nethome/axy148/R/lib64/R/library/ChipSeq/extdata/zhao_data.csv","/nethome/axy148/_visualization/sort_bam.txt","~/BamCompareZhao")'

useBamCompare2 <- function(input.sample.file,input.bam.file,output.dir)
{
  
  #sampleF <- read.table(input.sample.file,header = FALSE)
  #colnames(sampleF) <- c("ID","SampleName")
  
  #bam <- prepareBamFile(input.bam.file,file.pattern)
  
  re <- GetSampleInfo(input.sample.file, input.bam.file)
  
  cellInfo <- re$y
  # output.dir.name = dirname(input.sample.file)
  
  if (!dir.exists(output.dir))
  {
    dir.create(output.dir, recursive = TRUE)
  }
  
  temp3 = output.dir
  
  # cmd9 = 'ngs.plot.r -G' cmd10 = '-R' cmd11 = '-C' cmd12 = '-O' cmd13 = '-T'
  # cmd14 = '-L' cmd15 = '-RR' cmd16 = '-CD' cmd17= '-GO'
  
  
  
  cellInfo.run <- lapply(1:length(cellInfo), function(u,cellInfo, 
                                                      temp3)
  {
    
    x.name = cellInfo[[u]]$name
    
    es <- cellInfo[[u]]$es
    
    x.input <- es[es$Type_TF == "Input", ]$file.name
    
    x.sample <- es[es$Type_TF != "Input", ]
    
    #print(x.sample)
    #print(x.input)
    
    x.run <- apply(x.sample, 1, function(x, x.input, temp3)
    {
      
      y <- x
      
      ID <- y[1]
      Type_Cell <- y[2]
      Type_TF <- y[3]
      Cell_TF <- y[4]
      file.name <- y[5]
      xx <- file.name
      xx.name = paste(ID, gsub(" ", "-", Type_Cell), Type_TF, sep = "-")
      
      # ~/python/Python-2.7.11/python  ~/NGS_tools/deepTools/bin/bamCompare  --bamfile1 /scratch/projects/bbc/aiminy_project/DannyNewNgsPlot/2017-03-02-03_S11_R1.marked_sorted.bam --bamfile2 /scratch/projects/bbc/aiminy_project/DannyNewNgsPlot/2017-03-02-17_S1_R1.marked_sorted.bam --binSize 25 --ratio log2 -o ~/BamCompare/log2ratio_2017-03-02-03_S11_R1.marked_sorted.bam_vs_2017-03-02-17_S1_R1.marked_sorted.bam.bw 
      
      cmd1 <- paste("~/python/Python-2.7.11/python  ~/NGS_tools/deepTools/bin/bamCompare --bamfile1",xx,"--bamfile2",x.input,sep=" ")
      cmd2 <- "--binSize 25"
      cmd3 <- "--ratio log2 -o"
      
      cmd4 <- paste(cmd1,cmd2,cmd3,sep=" ") 
      
      cmd5 <- file.path(output.dir,paste0("log2ratio_",basename(as.character(xx)),"_vs_",basename(as.character(x.input)),".bw"))
      
      
      cmd6 <- paste(cmd4,cmd5,sep=" ") 
      
      cmd6
      
      #cat(cmd6, "\n")
      #cat("\n")
    }, x.input, temp3)
    
    x.run
    
  }, cellInfo,temp3)
  
  
  names(cellInfo.run) = unlist(lapply(cellInfo, function(u)
  {
    u$name
  }))
  
  system("module unload python/2.7.3")
  
  zzz <- unlist(cellInfo.run)
  
  lapply(1:length(zzz), function(u, zzz)
  {
    
    cat(as.character(zzz[u][[1]]), "\n")
    cat("\n")
    
    system(as.character(zzz[u][[1]]))
    
  }, zzz)
  # # dir.name=temp3 dir.name=reformatPath(dir.name)
  # 
  # file.name = file.path(temp3, dir(temp3, recursive = TRUE))
  # 
  # file.name.2 <- as.list(file.name)
  # 
  # 
  # # names(file.name.2) = unlist(lapply(file.name.2, function(u) { u$name }))
  # 
  # zzz <- unlist(file.name.2)
  # 
  # lapply(1:length(zzz), function(u, zzz)
  # {
  #   
  #   dir.name = dirname(zzz[u][[1]])
  #   file_name = file_path_sans_ext(basename(zzz[u][[1]]))
  #   
  #   cmd = paste("ngs.plot.r -G hg19 -R tss -C", zzz[u][[1]], "-O", file.path(dir.name, 
  #                                                                            paste0(file_name, ".tss")), "-T", file_name, "-L 4000 -RR 1 -CD 1 -GO total", 
  #               sep = " ")
  #   
  #   
  #   
  #   # system(as.character(zzz[u][[1]])) job.name = paste0('bamPlot.', u)
  #   # cmd.pegasus = usePegasus(job.option, Wall.time = '72:00',cores = 32,Memory
  #   # = 16000,span.ptile = 16,job.name) cmd2 = paste(cmd.pegasus,cmd,sep = ' ')
  #   
  #   cmd2 = cmd
  #   cat(cmd2, "\n")
  #   cat("\n")
  #   system(cmd2)
  # }, zzz)
  # 
  # 
  # re <- list(cellInforun = cellInfo.run, zzz = zzz)
  # 
  # # AnntationUsingChipSeeker(temp3, 'peaks.bed', temp3, DD = 5000)
  
  #return(re)
  
}

# bsub -P bbc -J "zhaoCJun" -o %J.zhaoCJun.log -e %J.zhaoCJun.err -W 72:00 -n 32 -q parallel -R 'rusage[mem= 16000 ] span[ptile= 16 ]' -u aimin.yan@med.miami.edu R -e 'library(ChipSeq);re <- ChipSeq:::ngs()'

ngs <- function(){
  
system("ngs.plot.r -G hg19 -R tss -C ~/zhao_config.txt -O zhao.tss -L 4000")

}

# bsub -P bbc -J "zhaoCJun" -o %J.zhaoCJun.log -e %J.zhaoCJun.err -W 72:00 -n 32 -q parallel -R 'rusage[mem= 16000 ] span[ptile= 16 ]' -u aimin.yan@med.miami.edu R -e 'library(ChipSeq);re <- ChipSeq:::ngs2()'

ngs2 <- function(config.file,distance.around.tss,output.file.dir){
  
  if (!dir.exists(output.file.dir))
  {
    dir.create(output.file.dir, recursive = TRUE)
  }
  
  cmd=paste("ngs.plot.r -G hg19 -R tss -C",config.file,"-O",output.file.dir,"-L",distance.around.tss,"-CO white:white:red", collapse = " ")
  
  system(cmd)
  
}

doAlignment<-function(input.fastq.dir,bowtie2Indeces,output.dir){
  re <- parserreadfiles(input.fastq.dir, "fastq")
  res <- re$input
  index <- system('echo $LSB_JOBINDEX',intern = TRUE)
  cmd0="bowtie2 -x"
  "/home/someusername/bowtie2Indeces/somegenome -U /home/someusername/myFastqFile.fastq -S /home/someusername/myFastqFile.sam"
  u <- as.integer(index)
  path_name = dirname(res[[u]])
  file_name = file_path_sans_ext(basename(res[[u]]))
  cmd1 <- paste(cmd0,bowtie2Indeces,"-U",res[[u]],"-S",file.path(output.dir,paste0(file_name,".sam")), sep = " ")
  cmd2 <- cmd1
  system(cmd2)
  cat(cmd2,"\n")
}

filterSam<-function(input.fastq.dir,bowtie2Indeces,output.dir){
  re <- parserreadfiles(input.fastq.dir, "fastq")
  res <- re$input
  index <- system('echo $LSB_JOBINDEX',intern = TRUE)
  cmd2='samtools view -Sh /home/someusername/myFastqFile.sam | grep -e "^@" -e "XM:i:[012][^0-9]" | grep -v "XS:i:" > /home/someusername/myFastqFile.sam.filtered.sam'
  
  cmd0="bowtie2 -x"
  "/home/someusername/bowtie2Indeces/somegenome -U /home/someusername/myFastqFile.fastq -S /home/someusername/myFastqFile.sam"
  u <- as.integer(index)
  path_name = dirname(res[[u]])
  file_name = file_path_sans_ext(basename(res[[u]]))
  cmd1 <- paste(cmd0,bowtie2Indeces,"-U",res[[u]],"-S",file.path(output.dir,paste0(file_name,".sam")), sep = " ")
  cmd2 <- cmd1
  system(cmd2)
  cat(cmd2,"\n")
}

convertSam2Bam<-function(input.fastq.dir,bowtie2Indeces,output.dir){
  re <- parserreadfiles(input.fastq.dir, "fastq")
  res <- re$input
  index <- system('echo $LSB_JOBINDEX',intern = TRUE)
  cmd3="samtools view -S -b /home/someusername/myFastqFile.sam.filtered.sam > /home/someusername/myFastqFile.sam.filtered.sam.bam"
  cmd0="bowtie2 -x"
  "/home/someusername/bowtie2Indeces/somegenome -U /home/someusername/myFastqFile.fastq -S /home/someusername/myFastqFile.sam"
  u <- as.integer(index)
  path_name = dirname(res[[u]])
  file_name = file_path_sans_ext(basename(res[[u]]))
  cmd1 <- paste(cmd0,bowtie2Indeces,"-U",res[[u]],"-S",file.path(output.dir,paste0(file_name,".sam")), sep = " ")
  cmd2 <- cmd1
  system(cmd2)
  cat(cmd2,"\n")
}

sortBam<-function(input.fastq.dir,bowtie2Indeces,output.dir){
  re <- parserreadfiles(input.fastq.dir, "fastq")
  res <- re$input
  index <- system('echo $LSB_JOBINDEX',intern = TRUE)
  cmd4="samtools sort /home/someusername/myFastqFile.sam.filtered.sam.bam /home/someusername/myFastqFile.sam.filtered.sam.bam.sorted"
  cmd0="bowtie2 -x"
  "/home/someusername/bowtie2Indeces/somegenome -U /home/someusername/myFastqFile.fastq -S /home/someusername/myFastqFile.sam"
  u <- as.integer(index)
  path_name = dirname(res[[u]])
  file_name = file_path_sans_ext(basename(res[[u]]))
  cmd1 <- paste(cmd0,bowtie2Indeces,"-U",res[[u]],"-S",file.path(output.dir,paste0(file_name,".sam")), sep = " ")
  cmd2 <- cmd1
  system(cmd2)
  cat(cmd2,"\n")
}

removePCRduplicates<-function(input.fastq.dir,bowtie2Indeces,output.dir){
  re <- parserreadfiles(input.fastq.dir, "fastq")
  res <- re$input
  index <- system('echo $LSB_JOBINDEX',intern = TRUE)
  
  cmd5="samtools rmdup -s /home/someusername/myFastqFile.sam.filtered.sam.bam.sorted.bam /home/someusername/myFastqFile.sam.filtered.sam.bam.sorted.bam.nodup.bam"
  cmd0="bowtie2 -x"
  "/home/someusername/bowtie2Indeces/somegenome -U /home/someusername/myFastqFile.fastq -S /home/someusername/myFastqFile.sam"
  u <- as.integer(index)
  path_name = dirname(res[[u]])
  file_name = file_path_sans_ext(basename(res[[u]]))
  cmd1 <- paste(cmd0,bowtie2Indeces,"-U",res[[u]],"-S",file.path(output.dir,paste0(file_name,".sam")), sep = " ")
  cmd2 <- cmd1
  system(cmd2)
  cat(cmd2,"\n")
}

indexBAM<-function(input.fastq.dir,bowtie2Indeces,output.dir){
  re <- parserreadfiles(input.fastq.dir, "fastq")
  res <- re$input
  index <- system('echo $LSB_JOBINDEX',intern = TRUE)
  
  cmd6="samtools index /home/someusername/myFastqFile.sam.filtered.sam.bam.sorted.bam.nodup.bam"
  
  cmd0="bowtie2 -x"
  "/home/someusername/bowtie2Indeces/somegenome -U /home/someusername/myFastqFile.fastq -S /home/someusername/myFastqFile.sam"
  u <- as.integer(index)
  path_name = dirname(res[[u]])
  file_name = file_path_sans_ext(basename(res[[u]]))
  cmd1 <- paste(cmd0,bowtie2Indeces,"-U",res[[u]],"-S",file.path(output.dir,paste0(file_name,".sam")), sep = " ")
  cmd2 <- cmd1
  system(cmd2)
  cat(cmd2,"\n")
}

bam2Bed<-function(input.fastq.dir,bowtie2Indeces,output.dir){
  re <- parserreadfiles(input.fastq.dir, "fastq")
  res <- re$input
  index <- system('echo $LSB_JOBINDEX',intern = TRUE)
  cmd7="bedtools bamtobed -i /home/someusername/myFastqFile.sam.filtered.sam.bam.sorted.bam.nodup.bam > /home/someusername/myFastqFile.sam.filtered.sam.bam.sorted.bam.nodup.bam.bed"
  cmd0="bowtie2 -x"
  "/home/someusername/bowtie2Indeces/somegenome -U /home/someusername/myFastqFile.fastq -S /home/someusername/myFastqFile.sam"
  u <- as.integer(index)
  path_name = dirname(res[[u]])
  file_name = file_path_sans_ext(basename(res[[u]]))
  cmd1 <- paste(cmd0,bowtie2Indeces,"-U",res[[u]],"-S",file.path(output.dir,paste0(file_name,".sam")), sep = " ")
  cmd2 <- cmd1
  system(cmd2)
  cat(cmd2,"\n")
}

getBamReady<-function(input.fastq.dir,bowtie2Indecesinput,sample.num,output.bam.dir){

Rfun1 <- 'library(ChipSeq);library(DoGs);re <- ChipSeq:::doAlignment('
Rinput <- paste0('\\"',input.fastq.dir,'\\",',
                 '\\"',bowtie2Indecesinput,'\\",',
                 '\\"',output.dir,'\\"')
Rfun2 <- ')'
Rfun <-paste0(Rfun1,Rinput,Rfun2)
alignment <- createBsubJobArrayRfun(Rfun,paste0("alignment[1-",sample.num,"]",NULL))
system(alignment)
}

getEffectOfRmDup<-function(){
  not.rm <- read.table("~/Dropbox (BBSR)/Aimin_project/Research/ChipSeq/Danny_ChipSeq_data_QC/qc_16.txt")
  rm <- read.table("~/Dropbox (BBSR)/Aimin_project/Research/ChipSeq/Danny_ChipSeq_data_QC_2/16_samples_QC.txt")
  res <- as.data.table(rbind(cbind(rep("not_rm",length(not.rm[,11])),not.rm[,11]),cbind(rep("rm",length(rm[,11])),rm[,11])))
  boxplot(as.numeric(res$V2)~res$V1)
}


#' sbatch -P bbc -J "zhaoCJun" -o %J.zhaoCJun.log -e %J.zhaoCJun.err -W 72:00 -n 32 -q parallel -R 'rusage[mem= 16000 ] span[ptile= 16 ]' -u aimin.yan@med.miami.edu R -e 'library(ChipSeq);re <- ChipSeq:::runIDR()'

#' sbatch -p hipri test-job.sh
#' sbatch -p hipri test-job.sh
#' squeue -u steven.miamiuniv


runIDR <- function()
{
cmd="idr --samples /home/steven.miamiuniv/Softwares/idr-2.0.2/tests/data/peak1 /home/steven.miamiuniv/Softwares/idr-2.0.2/tests/data/peak2 --plot testIDR"
system(cmd)

}

# bsub -P bbc -J "zhaoCJun" -o %J.zhaoCJun.log -e %J.zhaoCJun.err -W 72:00 -n 32 -q parallel -R 'rusage[mem= 16000 ] span[ptile= 16 ]' -u aimin.yan@med.miami.edu R -e 'library(ChipSeq);re <- ChipSeq:::peakCallAndAnnotationWithoutInput("/scratch/projects/bbc/aiminy_project/7_27_2017_bam/","/scratch/projects/bbc/aiminy_project/macs2_call","hs","macs2",0.0001)'
# 
peakCallAndAnnotationWithoutInput <- function(input.file.dir,output.file.dir,genome = "hs",peakcaller = c("macs14", "macs2"), peakPvalue) {
  
  genome <- match.arg(genome)
  
  cmd10 <- paste("-f BAM", "-g", genome, "-n", sep = " ")
  
  switch(peakcaller, macs2 = {
    
    PATH1 = Sys.getenv("PATH")
    
    macs2_Lib = file.path("/nethome/axy148/NGS_tools/MACS/bin/")
    
    Sys.setenv(PATH = paste0(macs2_Lib, ":", PATH1))
    
    cmd1 <- Sys.which("macs2")[[1]]
    
    cat(cmd1, "\n")
    
    cmd9 = paste(cmd1, "callpeak -t", sep = " ")
    
    cmd11 <- paste("-p", peakPvalue, sep = " ")
    
  }, {
    cmd9 = "macs14 -t "
    
    cmd11 <- paste("-m 6,18 --bw=200 -p", peakPvalue, sep = " ")
    
  })
  
  re<-ParserReadFiles(input.file.dir,"bam")
  
  file.name.2<-re$input
  #output.dir.name=re$output
  
  temp3=file.path(output.file.dir,"PeakCall")
  
  if(!dir.exists(temp3)){dir.create(temp3,recursive = TRUE)}
  
  re.out<-file.name.2
  
  #cmd9="macs14 -t "
  #cmd10="-f BAM -g hs -n "
  #cmd11=" -m 6,18 --bw=200 -p 0.00001"
  
  lapply(1:length(re.out),function(u,re.out,temp3){
    
    x=re.out[[u]]
    x_name=names(re.out)[u]
    
    cmd12=paste(cmd9,x,cmd10,file.path(temp3,paste0(x_name,"_hs_1.00e-05_macs")),cmd11,sep=" ")
    
    print(cmd12)
    system(cmd12, intern = TRUE, ignore.stderr = TRUE)
    
    #re=read.table(u,header=FALSE)
    #  re<-as.character(re[,1])
    #  #colnames(re)=c("Count","GeneName")
    #  re
  },re.out,temp3)
  
  #AnnotatePeak2(paste0(temp3,"/"),"*macs142_peaks.bed",7,paste0(output.dir.name,"PeakAnnotation_at_",temp2),genome="Hs")
  
  #AnnotatePeak3(paste0(temp3,"/"),paste0(output.dir.name,"_PeakAnnotation"),
  #              genome="Hs")
  
  #BamFileSortIndexVisualization2(re,genome)
  
}

#input <- "~/nodup.bam.txt"
#output <- "~/"
#
# R -e 'library(ChipSeq);library(ChipSeq);ChipSeq:::parseToSampleInfo("/projects/scratch/bbc/Project/Danny_chip3/Filtered_bam","*.bam$","~/Danny_chip3","sample_infor_Danny_chip3.txt","Danny_chip3")'
# 
parseToSampleInfo <- function(input.dir,input.pattern,output.dir,output.file,label){
  
  file.1 <- list.files(input.dir,pattern=input.pattern, all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
  s <- file.1
  ss <- basename(as.character(s))
  sss <- tools::file_path_sans_ext(ss)
  ssss <- tools::file_path_sans_ext(sss)
  
  sm.info <- cbind.data.frame(sss,rep(label,length(sss)),ssss,rep("Link",length(sss)))
  
  colnames(sm.info) <- c("ID","Type_Cell","Type_TF","Link")
  
  if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}
  
  write.table(sm.info,file=file.path(output.dir,output.file),row.names = FALSE,quote=FALSE,sep=",")
  
}

# ChipSeq:::getAb2inputPair("/Volumes/Bioinformatics$/Aimin_project/Danny_chip3.txt")
# 
getAb2inputPair <- function(input.infor.file)
{
  
  f <- read.table(input.infor.file,header = FALSE)
  
  colnames(f)=c("Ab","Input")
  
  print(f)
  
  return(f)
}
 
#  ChipSeq:::getBwUseBamCompare("/Volumes/Bioinformatics$/Aimin_project/Danny_chip3.txt","~")
#  
getBwUseBamCompare <- function(input.infor.file,bl_file,output.dir)
{
  x.sample <- getAb2inputPair(input.infor.file)
    
  if (!dir.exists(output.dir))
  {
    dir.create(output.dir, recursive = TRUE)
  }
  
    x.run <- apply(x.sample, 1, function(x,output.dir)
    {
      
      print(t(x))
      
      x <- as.data.frame(t(x))
      xx <- x$Ab
      xx.input <- x$Input
      
      #xx.name = paste(ID, gsub(" ", "-", Type_Cell), Type_TF, sep = "-")
      
      cmd1 <- paste("bamCompare --bamfile1",xx,"--bamfile2",xx.input,sep=" ")
      cmd2 <- "--binSize 25"
      #cmd3 <- "--ratio log2 -o"
      cmd3 <- "--ratio subtract --normalizeTo1x 2451960000 --ignoreForNormalization chrX -o"
      cmd3 <- "--ratio subtract --normalizeTo1x 2451960000 --ignoreForNormalization chrX --blackListFileName"
      cmd3 <- paste(cmd3,bl_file,"-o",sep=" ")
      cmd4 <- paste(cmd1,cmd2,cmd3,sep=" ") 
      
      cmd5 <- file.path(output.dir,paste0("log2ratio_",basename(as.character(xx)),"_vs_",basename(as.character(xx.input)),".bw"))
      
      cmd6 <- paste(cmd4,cmd5,sep=" ") 
      
      cmd6
      
      print(cmd6)
      
      system(cmd6)
      
    }, output.dir)
    

}

# R -e 'library(ChipSeq)::peakCall2("/Volumes/Bioinformatics$/Aimin_project/Danny_chip3.txt","hs","~/","macs2",0.00001)'
# 
peakCall2 <- function(input.infor.file,genome = c("Hs","hs", "HS", "hS"),output.dir,peakcaller = c("macs14","macs2"), peakPvalue)
{
  x.sample <- getAb2inputPair(input.infor.file)
  
  if (!dir.exists(output.dir))
  {
    dir.create(output.dir, recursive = TRUE)
  }
  
  peakcaller <- match.arg(peakcaller)

  genome <- match.arg(genome)

  genome <- tolower(genome)
  
  cmd10 <- paste("-f BAM", "-g", genome, "-n", sep = " ")
  
  switch(peakcaller, macs2 = {
    
    cmd1 <- "macs2"
    
    cat(cmd1, "\n")
    
    cmd9 = paste(cmd1, "callpeak -t", sep = " ")
    
    cmd11 <- paste("-p", peakPvalue, sep = " ")
    
  }, {
    cmd9 = "macs14 -t "
    
    cmd11 <- paste("-p", peakPvalue, sep = " ")
    
  })
  
  macs.run <- apply(x.sample,1,function(u, output.dir)
  {
    x <- as.data.frame(t(u))
    xx <- x$Ab
    xx.input <- x$Input
    
    x.name <- paste(peakcaller,basename(as.character(xx)),basename(as.character(xx.input)),genome,peakPvalue,sep="_")
    
    cmd12 = paste(cmd9, xx, "-c", xx.input, cmd10, file.path(output.dir, x.name), cmd11, sep = " ")
      
    print(cmd12)
      
    system(cmd12)
    
    },output.dir)

  #AnntationUsingChipSeeker(output.dir, "peaks.bed", temp3, DD = 5000)
  
#  return(re)
  
}

useBamCoverage <- function(input.sample.file,input.bam.file,output.dir)
{
  
  re <- GetSampleInfo(input.sample.file, input.bam.file)
  
  cellInfo <- re$y
  
  # output.dir.name = dirname(input.sample.file)
  
  if (!dir.exists(output.dir))
  {
    dir.create(output.dir, recursive = TRUE)
  }
  
  temp3 = output.dir
  
  # cmd9 = 'ngs.plot.r -G' cmd10 = '-R' cmd11 = '-C' cmd12 = '-O' cmd13 = '-T'
  # cmd14 = '-L' cmd15 = '-RR' cmd16 = '-CD' cmd17= '-GO'
  
  
  
  cellInfo.run <- lapply(1:length(cellInfo), function(u,cellInfo, 
                                                      temp3)
  {
    
    x.name = cellInfo[[u]]$name
    
    es <- cellInfo[[u]]$es
    
    x.input <- es[es$Type_TF == "Input", ]$file.name
    
    x.sample <- es[es$Type_TF != "Input", ]
    
    #print(x.sample)
    #print(x.input)
    
    x.run <- apply(x.sample, 1, function(x, x.input, temp3)
    {
      
      y <- x
      
      ID <- y[1]
      Type_Cell <- y[2]
      Type_TF <- y[3]
      Cell_TF <- y[4]
      file.name <- y[5]
      xx <- file.name
      xx.name = paste(ID, gsub(" ", "-", Type_Cell), Type_TF, sep = "-")
      
      # ~/python/Python-2.7.11/python  ~/NGS_tools/deepTools/bin/bamCompare  --bamfile1 /scratch/projects/bbc/aiminy_project/DannyNewNgsPlot/2017-03-02-03_S11_R1.marked_sorted.bam --bamfile2 /scratch/projects/bbc/aiminy_project/DannyNewNgsPlot/2017-03-02-17_S1_R1.marked_sorted.bam --binSize 25 --ratio log2 -o ~/BamCompare/log2ratio_2017-03-02-03_S11_R1.marked_sorted.bam_vs_2017-03-02-17_S1_R1.marked_sorted.bam.bw 
      
      #bamCoverage -b reads.bam -o coverage.bw
      
      #bamCoverage --bam a.bam -o a.SeqDepthNorm.bw \
      #--binSize 10
      #--normalizeTo1x 2150570000
      #--ignoreForNormalization chrX
      #--extendReads
      
      cmd1 <- paste("~/python/Python-2.7.11/python  ~/NGS_tools/deepTools/bin/bamCoverage --bam",xx,sep=" ")
      cmd2 <- "--binSize 25 --normalizeTo1x 2451960000 --ignoreForNormalization chrX"
      cmd3 <- "-o"
      
      cmd4 <- paste(cmd1,cmd2,cmd3,sep=" ") 
      
      cmd5 <- file.path(output.dir,paste0("To1x_",basename(as.character(xx)),"_vs_","no_input",".bw"))
      
      
      cmd6 <- paste(cmd4,cmd5,sep=" ") 
      
      cmd6
      
      #cat(cmd6, "\n")
      #cat("\n")
    }, x.input, temp3)
    
    x.run
    
  }, cellInfo,temp3)
  
  
  names(cellInfo.run) = unlist(lapply(cellInfo, function(u)
  {
    u$name
  }))
  
  zzz <- unlist(cellInfo.run)
  
  lapply(1:length(zzz), function(u, zzz)
  {
    
    cat(as.character(zzz[u][[1]]), "\n")
    cat("\n")
    
    system(as.character(zzz[u][[1]]))
    
  }, zzz)
  # # dir.name=temp3 dir.name=reformatPath(dir.name)
  # 
  # file.name = file.path(temp3, dir(temp3, recursive = TRUE))
  # 
  # file.name.2 <- as.list(file.name)
  # 
  # 
  # # names(file.name.2) = unlist(lapply(file.name.2, function(u) { u$name }))
  # 
  # zzz <- unlist(file.name.2)
  # 
  # lapply(1:length(zzz), function(u, zzz)
  # {
  #   
  #   dir.name = dirname(zzz[u][[1]])
  #   file_name = file_path_sans_ext(basename(zzz[u][[1]]))
  #   
  #   cmd = paste("ngs.plot.r -G hg19 -R tss -C", zzz[u][[1]], "-O", file.path(dir.name, 
  #                                                                            paste0(file_name, ".tss")), "-T", file_name, "-L 4000 -RR 1 -CD 1 -GO total", 
  #               sep = " ")
  #   
  #   
  #   
  #   # system(as.character(zzz[u][[1]])) job.name = paste0('bamPlot.', u)
  #   # cmd.pegasus = usePegasus(job.option, Wall.time = '72:00',cores = 32,Memory
  #   # = 16000,span.ptile = 16,job.name) cmd2 = paste(cmd.pegasus,cmd,sep = ' ')
  #   
  #   cmd2 = cmd
  #   cat(cmd2, "\n")
  #   cat("\n")
  #   system(cmd2)
  # }, zzz)
  # 
  # 
  # re <- list(cellInforun = cellInfo.run, zzz = zzz)
  # 
  # # AnntationUsingChipSeeker(temp3, 'peaks.bed', temp3, DD = 5000)
  
  #return(re)
  
}

getBwUseBamCoverage2 <- function(input.infor.file,bl_file,output.dir)
{
  x.sample <- getAb2inputPair(input.infor.file)
  
  if (!dir.exists(output.dir))
  {
    dir.create(output.dir, recursive = TRUE)
  }
  
  x.run <- apply(x.sample, 1, function(x,output.dir)
  {
    
    print(t(x))
    
    x <- as.data.frame(t(x))
    xx <- x$Ab
    xx.input <- x$Input
    
    #xx.name = paste(ID, gsub(" ", "-", Type_Cell), Type_TF, sep = "-")
    
    cmd1 <- paste("bamCoverage --bam",xx,sep=" ")
    cmd2 <- "--binSize 25 --normalizeTo1x 2451960000 --ignoreForNormalization chrX --blackListFileName"
    cmd2 <- paste(cmd2,bl_file,sep=" ")
    cmd3 <- "-o"
    
    
    
    #cmd1 <- paste("bamCompare --bamfile1",xx,"--bamfile2",xx.input,sep=" ")
    #cmd2 <- "--binSize 25"
    #cmd3 <- "--ratio log2 -o"
    
    cmd4 <- paste(cmd1,cmd2,cmd3,sep=" ") 
    
    cmd5 <- file.path(output.dir,paste0("Un_norma",basename(as.character(xx)),".bw"))
    
    cmd6 <- paste(cmd4,cmd5,sep=" ") 
    
    cmd6
    
    print(cmd6)
    
    system(cmd6)
    
  }, output.dir)
  
  
}

getSummitSequence<-function(dir.name,input.file.pattern,genome,out.dir.name){
  
  if (!dir.exists(out.dir.name))
  {
    dir.create(out.dir.name, recursive = TRUE)
  }
  
  file.1 <- list.files(dir.name,pattern=input.file.pattern,all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
  
  print(file.1)
  
  re.out<-lapply(file.1,function(u){
    
    #print(names(u))
    uu=u
    re=read.table(uu,header=TRUE)
    
    print(head(re))
    
    colnames(re)[c(7,9)]=c("-10*LOG10(pvalue)","FDR(%)")
    temp<-re
    temp$start<-temp$start-1
    temp$end<-temp$end-1
    temp$summit<-temp$abs_summit-1
    temp
    
    summitPeak<-temp$start+temp$summit
    
    temp2<-temp
    
    temp2$start<-summitPeak-49
    temp2$end<-summitPeak+50
    temp2$summit<-summitPeak
    temp2$length<-temp2$end-temp2$start+1
    temp2
    rownames(temp2)<-paste0("MACS_peak_",rownames(temp2))
    temp3<-toGRanges(temp2)
    names(temp3)<-rownames(temp2)
    
    
    if(genome == "mm10")
    {
    dd.GRCm39.mm10<-toGRanges(EnsDb.Mmusculus.v75)
    genome(temp3)<-genome(dd.GRCm39.mm10)
    seqlevels(dd.GRCm39.mm10,force=TRUE) <- seqlevels(temp3)
    seqinfo(temp3)<-seqinfo(dd.GRCm39.mm10)
    temp3.trimmed<-trim(temp3, use.names=TRUE)
    genome.mm10<-getBSgenome("BSgenome.Mmusculus.UCSC.mm10")
    genome(temp3.trimmed)<-"mm10"
    }
    
    if(genome == "hg19"){
      
      dd.GRCm39.mm10<-toGRanges(EnsDb.Hsapiens.v75)
      genome(temp3)<-genome(dd.GRCm39.mm10)
      seqlevels(dd.GRCm39.mm10,force=TRUE) <- seqlevels(temp3)
      seqinfo(temp3)<-seqinfo(dd.GRCm39.mm10)
      temp3.trimmed<-trim(temp3, use.names=TRUE)
      genome.mm10<-getBSgenome("BSgenome.Hsapiens.UCSC.hg19")
      genome(temp3.trimmed)<-"hg19"
      
    }
    
    
    seq.temp3<-getAllPeakSequence(temp3.trimmed,genome=genome.mm10)
    
    write2FASTA(seq.temp3, file.path(out.dir.name,paste0(basename(u),".fa")))
    
    #re2<-list(originalPeak=re,down1Peak=temp,aroundSummit100Peak=temp2,GR=temp3)
    #re2
  })
  
  #sample.name<-sapply(strsplit(names(file.name.2),split="_peaks_"),"[[",1)
  
  #names(re.out)<-sample.name
  
  
  #return(re.out)
}

#' input.bed.file <- '/Users/aiminyan/Aimin/DropboxUmass/NADfinder/BedFiles/"Peric-Hupkes 2010 MEF LADs mm10.bed"'
#' map.chain.file <- "~/Aimin/DropboxUmass/Aimin/Project/ReferenceGenome/mm10ToMm9.over.chain"
#' hmm.file <- "/Volumes/Aimin4TB-2/Aimin_Project/Ubuntu_local/chromatin_states_chromHMM_mm9/spleen_cStates_HMM.bed"
#' useChromHMM(input.bed.file,map.chain.file,hmm.file,output.file.dir)

useChromHMM <- function(input.bed.file,map.chain.file,hmm.file,output.file.dir) {
  
  x_name <- gsub(" ","-",tools::file_path_sans_ext(basename(input.bed.file)))
  
  output.dir <- file.path(output.file.dir,"histone_modification")
  
  if (!dir.exists(output.dir)){dir.create(output.dir, recursive = TRUE)}
  
  temp1 <- tempfile()
  temp2 <- tempfile()
  temp3 <- tempfile()
  
  cmd0 = paste("liftOver",input.bed.file,map.chain.file,temp1,temp2)
  system(cmd0)

  peak.all.in.mm9 <- read.table(file = temp1)
  peak.all.in.mm9 <- peak.all.in.mm9[,1:3]
  
  colnames(peak.all.in.mm9) <- c("chr","start","end")
  
  cmd1 = paste("bedtools intersect -a",hmm.file,"-b",temp1,">",temp3)
  system(cmd1)  

  peakAll.anno <- read.table(temp3)
  
  head(peakAll.anno)
  
  ID <- paste(as.character(peakAll.anno[,1]),peakAll.anno[,2],peakAll.anno[,3],sep = "_")
  
  peakAll.anno.2 <- cbind.data.frame(peakAll.anno,ID=ID)
  
  ID2 <- paste(as.character(peak.all.in.mm9[,1]),peak.all.in.mm9[,2],peak.all.in.mm9[,3],sep = "_")
  
  peakAll.anno.3 <- cbind.data.frame(peak.all.in.mm9,ID=ID2)
  
  peakAll.anno.4 <- merge.data.frame(peakAll.anno.2,peakAll.anno.3,by="ID")
  
  getState <- function(g){
    
    y <- lapply(g, function(u){
      
      if(!is.na(str_locate(u,"_")[1])){
        x <- str_sub(u,1,str_locate(u,"_")[1]-1)
      }else
      {
        x <- u
      }
      x
    })
    
    gg <- unlist(y)
    
    ggg <- gg[!is.na(gg)]
    ggg
    
  }
  
  getAnno <- function(g){
    
    y <- lapply(g, function(u){
      
      if(!is.na(str_locate(u,"_")[1])){
        x <- str_sub(u,str_locate(u,"_")[1]+1,str_length(u))
      }else
      {
        x <- u
      }
      x
    })
    
    gg <- unlist(y)
    
    ggg <- gg[!is.na(gg)]
    ggg
    
  }
  
  state <- getState(peakAll.anno.4$V4)
  anno <- getAnno(peakAll.anno.4$V4)
  
  anno.6 <- anno
  
  anno.6[which(anno.6 %in% c("Active_Promoter","Poised_Promoter"))] <- "Promoter"
  
  anno.6[which(anno.6 %in% c("Strong_Enhancer","Poised_Enhancer"))] <- "Enhancer"
  
  anno.6[which(anno.6 %in% c("Txn_Transition","Txn_Elongation","Weak_Txn"))] <- "Transcribed"
  
  peakAll.anno.5 <- cbind.data.frame(peakAll.anno.4,state=state,anno=anno,anno.6 = anno.6)
  
  mytable <- table(peakAll.anno.5$anno.6)
  
  lbls <- paste(names(mytable),mytable, sep="-")
  pie(mytable, labels = lbls,main="Annotation")
  
  mytable2 <- as.data.frame(mytable)
  colnames(mytable2)= c("Function_Annotation","counts")
  
  table_lables <- mytable2 %>%
    mutate(Function_Annotation=factor(Function_Annotation,levels=Function_Annotation[length(Function_Annotation):1]),
           cumulative=cumsum(counts),
           midpoint= cumulative-(counts/2),
           labels=paste0(round((counts/sum(counts))*100,2),"%"," (",counts,") "))
  
  sp1 <- ggplot(table_lables,aes(x="",y=counts,fill=Function_Annotation))+
    geom_bar(width = 1,stat="identity") +
    coord_polar(theta="y",start = 0,direction = 1) +
    scale_fill_manual(values = c("Lightblue","#AD7366","Lightgreen","Orange","Coral","Yellow"))+
    labs(x="",y="",title=x_name,fill="Function_Annotation")+
    geom_text(aes(x=1.2,y=midpoint,label=labels),color="black",fontface="bold",size=3.3) +
    theme(plot.title=element_text(hjust=0.5),
          legend.title = element_text(hjust = 0.5,face="bold",size=10))
  
  spp <- list(sp1=sp1)
  
  multi.page <- ggarrange(plotlist=spp,nrow = 1, ncol = 1)
  
  if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}
  ggexport(multi.page, filename = file.path(output.dir,paste0(x_name,"_annotation.pdf")))
}

# input.file <- "~/Downloads/files.txt"
downLoadHMMChipSeqData <- function(input.file) {
  cmd ="xargs -n 1 curl -O -L < input.file"
  system(cmd)
}

#  input.file.dir <- "~/Aimin/ProjectAtCluster/umw_ingolf_bach/Aimin/MergedFq"
#  input.file.pattern <- "*fastq$"
#  output.file.dir <- "~/Aimin/ProjectAtCluster/umw_ingolf_bach/Aimin/MergedFq"
#  output.file.name <- "fq4alignment.txt"
#  
#  readFqFiles(input.file.dir,input.file.pattern,output.file.dir,output.file.name)

readFqFiles <- function(input.file.dir,input.file.pattern,output.file.dir,output.file.name) {
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  file.name.4 <- list.files(input.file.dir,pattern=input.file.pattern,all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
  
  x <- basename(file.name.4)
  xx <- unique(str_sub(x,1,str_locate(x,"_R")[,1]-1))
  
  dfL <- lapply(xx, function(u,file.name.4){
    
    y <- basename(file.name.4)
    yy <- str_sub(x,1,str_locate(x,"_R")[,1]-1)
    
    yyy <- data.frame(sampeName=u,fq1=file.name.4[c(which(u==yy))][1],fq2=file.name.4[c(which(u==yy))][2])
    yyy
  
  },file.name.4)
  
  DF <- do.call(rbind.data.frame,dfL)
  
  write.table(DF,file = file.path(output.file.dir,output.file.name),append = FALSE, quote = F, sep = ",",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = F)
}

plotVenn<-function(venn_cnt, vennx, otherCounts=NULL, 
                   cat.cex = 1, cat.col = "black", 
                   cat.fontface = "plain", cat.fontfamily = "serif", ...){
  op=par(mar=c(0,0,0,0))
  on.exit(par(op))
  plot.new()
  venngrid <- venn.diagram(x=vennx, filename=NULL, cat.cex = cat.cex,
                           cat.col = cat.col, cat.fontface = cat.fontface, 
                           cat.fontfamily = cat.fontfamily, ...)
  unlink(dir(pattern="^VennDiagram[0-9_\\-]+.log$")) ## delete the log file
  if(grepl("^count\\.", colnames(venn_cnt)[ncol(venn_cnt)]) && 
     connectedPeaks=="keepAll"){
    n <- which(colnames(venn_cnt)=="Counts")-1
    venn_cnt.Counts.gt.1 <- rowSums(venn_cnt[ ,1:n]) > 1
    counts <- venn_cnt[-1,"Counts"]
    counts <- counts[counts!=0]
    if(!any(duplicated(counts))){## there is no duplicated counts number
      w <- 1
      for(i in 1:length(venngrid)){
        if(inherits(venngrid[[i]], what="text") && w<nrow(venn_cnt)){
          w <- w+1
          cnt <- as.numeric(venngrid[[i]]$label)
          if(!is.na(cnt) && cnt!=0){
            j <- which(venn_cnt[,n+1]==cnt)
            if(venn_cnt.Counts.gt.1[j]){
              labels <- venn_cnt[j, -(1:(n+1))]
              if(nrow(venn_cnt)==4 && sum(venn_cnt[,4]) < 
                 sum(venn_cnt[,5])){
                labels <- rev(labels)
              }
              labels <- paste(as.character(labels[labels>0]),
                              collapse="/", sep="")
              if(labels!="") venngrid[[i]]$label <- 
                paste(venngrid[[i]]$label, "(", 
                      labels, ")", sep="")
            }
          }
        }
      }
    }else{##otherwise guess the order of the number, 
      ##Bugs, the map order maybe wrong
      j <- 1
      map <- switch(as.character(nrow(venn_cnt)),
                    '4'=c(3, 4, 2),
                    '8'=c(5, 7, 3, 6, 8, 4, 2),
                    '16'=c(3, 4, 2, 11, 12, 16, 8, 6, 9, 10, 
                           14, 15, 7, 5, 13),
                    '32'=c(17, 9, 5, 3, 2, 6, 18, 19, 25, 10, 13, 21, 7, 
                           11, 4, 8, 22, 20, 27, 26, 30, 29, 23, 
                           15, 12, 16, 24, 28, 30, 31, 32))
      for(i in 1:length(venngrid)){
        if(inherits(venngrid[[i]], what="text") && j<length(map)){
          while(j < length(map) && as.numeric(venngrid[[i]]$label) != 
                venn_cnt[map[j],n+1]) j <- j+1
          if(as.numeric(venngrid[[i]]$label) == venn_cnt[map[j],
                                                         n+1] && 
             j <= length(map)){
            if(venn_cnt.Counts.gt.1[map[j]]){
              labels <- venn_cnt[map[j], -(1:(n+1))]
              labels <- paste(as.character(labels[labels>0]), 
                              collapse="/", sep="")
              if(labels!="") 
                venngrid[[i]]$label <- 
                paste(venngrid[[i]]$label, "(", 
                      labels, ")", sep="")
            }
            j <- j+1
          }
        }
      }
    }
  }
  if(!is.null(otherCounts)){
    tmp <- textGrob(label=otherCounts, x=0.9, y=0.1, 
                    gp=gpar(col = cat.col, cex = cat.cex, 
                            fontface = cat.fontface, 
                            fontfamily = cat.fontfamily))
    venngrid <- gList(venngrid, tmp)
  }
  grid.draw(venngrid)
}

testvenneuler <- function(){
A <- c("gene1", "gene2", "gene3", "gene5", "gene12", "", "")
B <- c("gene1", "gene2", "gene6", "gene7", "", "", "")
C <- c("gene2", "gene6", "gene7", "gene8", "gene9", "gene13", "gene14")
D <- c("gene7", "gene8", "gene9", "gene10", "gene11", "gene12", "")
dat <- data.frame(A,B,C,D)

vennfun <- function(x) { 
  x$id <- seq(1, nrow(x))  #add a column of numbers (required for melt)
  xm <- melt(x, id.vars="id", na.rm=TRUE)  #melt table into two columns (value & variable)
  xc <- dcast(xm, value~variable, fun.aggregate=length)  #remove NA's, list presence/absence of each value for each variable (1 or 0)
  rownames(xc) <- xc$value  #value column = rownames (required for Venneuler)
  xc$value <- NULL  #remove redundent value column
  xc  #output the new dataframe
}

library(reshape2)
install.packages("rJava")
library(rJava)
install.packages("venneuler")
library(venneuler)

VennDat <- vennfun(dat)
genes.venn <- venneuler(VennDat)
plot(genes.venn)

}

testColofulVenn <- function(){

  y <- c(6,856,297)
  names(y) <- c("01","10","11")
  labels <- c("TF1","TF2")
  plot.new()
  colorfulVennPlot::plotVenn(y, labels, Colors=rainbow(7))
}

orderPeakAndOutPut <- function(peakAll,output.file.dir,output.file.name,outHeader= FALSE) {
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  chrOrder<-c(paste("chr",1:19,sep=""),"chrX","chrY")
  
  peakAll$seqnames <- factor(peakAll$seqnames, levels=chrOrder)
  
  peakAll2 <- peakAll[order(peakAll$seqnames,peakAll$start),]
  
  if(outHeader == TRUE){
    write.table(peakAll2,file = file.path(output.file.dir,paste0(output.file.name,".txt")),append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T)
  }else{
    write.table(peakAll2,file = file.path(output.file.dir,paste0(output.file.name,".bed")),append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = F)
  }
}

# 
getPeaksInUnique <- function(ol.ciLAD.XL.MEF_LAD) {
  
  peaks.name <- names(ol.ciLAD.XL.MEF_LAD$uniquePeaks)
  uniquePeaks.name <- unique(str_sub(peaks.name,1,str_locate(peaks.name,"__")[,1]-1))
  
  Up <- lapply(uniquePeaks.name, function(u,ol.ciLAD.XL.MEF_LAD){
    
    unique.peaks <- ol.ciLAD.XL.MEF_LAD$uniquePeaks[grep(u,names(ol.ciLAD.XL.MEF_LAD$uniquePeaks)),]
    unique.peaks
  },ol.ciLAD.XL.MEF_LAD)
  names(Up) <- uniquePeaks.name
  
  cat("peak names in unique regions are: ",names(Up),"\n")
  
  cat('The peak you select:', '\n')
  
  m<-scan("",quiet=TRUE)
  
  unique.part.l <- unlist(lapply(1:length(m), function(u,m,Up){
    x= m[u]
    unique.part = Up[[x]]
    unique.part
  },m,Up))

  unique.part.name.l <- unlist(lapply(1:length(m), function(u,m,Up){
    x= m[u]
    unique.part.name = names(Up)[x]
    unique.part.name
  },m,Up))
  
  
  names(unique.part.l) <- unique.part.name.l
  
  unique.part.l
  
}

# XL.3.subsets <- getUniquePeaks(ol.ciLAD.XL.MEF_LAD)

# Given a list of peaks, to get enriched GO and pathway

# ol.ciLAD.XL.MEF_LAD <- findOverlapsOfPeaks(re.out[c(1,4,3)],connectedPeaks = "keepAll")
#
# 
# XL.3.subsets <- getUniquePeaks(ol.ciLAD.XL.MEF_LAD)
# output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_9_17_2018"
# genome = "Mm"
# 
# null <- peaksToEnrichedGO(XL.3.subsets,output.file.dir,genome)

# output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_10_22_2018/XL"
# genome = "Mm"
# null <- peaksToEnrichedGO(XL.3.subsets,output.file.dir,genome)
# output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_10_22_2018/XL/Bed"
# outGrl(XL.3.subsets,output.file.dir)

# output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_10_22_2018/nonXL"
# genome = "Mm"
# null <- peaksToEnrichedGO(nonXL.3.subsets,output.file.dir,genome)
# output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_10_22_2018/nonXL/Bed"
# outGrl(nonXL.3.subsets,output.file.dir)

# ol.ciLAD.non_XL.MEF_LAD <- findOverlapsOfPeaks(re.out[c(1,2,3)])

# nonXL.3.subsets <- getUniquePeaks(ol.ciLAD.non_XL.MEF_LAD)
# names(nonXL.3.subsets) <- gsub("XL","nonXL",names(nonXL.3.subsets))
# 
# output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_9_17_2018/nonXL"
# null <- peaksToEnrichedGO(nonXL.3.subsets,output.file.dir,genome)
 

# XL.3.subsets <- getUniquePeaks(ol.ciLAD.XL.MEF_LAD)
# output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_9_18_2018/XL"
# genome = "Mm"
# null <- peaksToEnrichedGO(XL.3.subsets,output.file.dir,genome)

# nonXL.3.subsets <- getUniquePeaks(ol.ciLAD.non_XL.MEF_LAD)
# output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_9_18_2018/nonXL"
# null <- peaksToEnrichedGO(nonXL.3.subsets,output.file.dir,genome)

peaksToEnrichedGO <- function(re.out,output.file.dir,genome){

  
 # ol.ciLAD.XL.MEF_LAD <- findOverlapsOfPeaks(re.out[c(1,4,3)],connectedPeaks = "keepAll")
  
  
  
# genome = "Mm"
null <- lapply(1:length(re.out),function(u,re.out,output.file.dir,genome){
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  #u <- 1
  #re.out <- XL.3.subsets
  
  x=re.out[[u]]
  x_name=names(re.out)[u]
  
  if(genome=="Mm"){
    
    #annoData <- toGRanges(EnsDb.Mmusculus.v75, feature="gene")
    
    dd.GRCm39.mm10<-toGRanges(EnsDb.Mmusculus.v75)
    
    overlaps.trimmed<-trim(x,use.names=TRUE)
    overlaps.anno<-annoPeaks(overlaps.trimmed,dd.GRCm39.mm10)
    
    #overlaps.anno<-annoPeaks(overlaps.trimmed,annoData)
    
    write.table(overlaps.anno,file=file.path(output.file.dir,paste0(x_name,"_annotation.txt")),row.names = FALSE,quote=FALSE,sep="\t")
    
     #annotatedPeak1 <- annotatePeakInBatch(overlaps.trimmed, AnnotationData=annoData)
    #  overlaps.anno.with.entrez.id <- addGeneIDs(annotatedPeak1,"org.Mm.eg.db",IDs2Add = "symbol")
    #  write.csv(as.data.frame(unname(overlaps.anno.with.entrez.id)), file.path(output.file.dir,paste0(x_name,"_other_anno.csv")))
    #  pdf(file.path(output.file.dir,paste0(x_name,"_annotation_pie_plot.pdf")))  
    #  pie1(table(overlaps.anno.with.entrez.id$insideFeature))
    #  dev.off()
    #  write.table(annotatedPeak1,file=file.path(output.file.dir,paste0(x_name,"_annotationInBatch.txt")),row.names = FALSE,quote=FALSE,sep="\t")
    
    getGoAndPath <- function(overlaps.anno.with.entrez.id,output.file.dir,x_name) {
      
      overlaps.anno.with.entrez.id <- overlaps.anno
      over <- getEnrichedGO(overlaps.anno.with.entrez.id, orgAnn="org.Mm.eg.db",
                            maxP=0.5, minGOterm=10,
                            multiAdjMethod="BH", condense=TRUE)
      
      
      path <- getEnrichedPATH(overlaps.anno.with.entrez.id, "org.Mm.eg.db", "reactome.db", maxP=.05)
      
      write.table(path,file=file.path(output.file.dir,paste0(x_name,"_path.txt")),row.names = FALSE,quote=FALSE,sep="\t")
      
      convert2geneSymbol <- function(over1) {
        
        geneSymbol <- lapply(over1$EntrezID, function(u){
          
          x <- annotate::getSYMBOL(unlist(strsplit(as.character(u),";")), data='org.Mm.eg') 
          x
          
        })
        
        list_to_df <- function(list_for_df)
        {
          list_for_df <- as.list(list_for_df)
          
          nm <- names(list_for_df)
          if (is.null(nm)) 
            nm <- seq_along(list_for_df)
          
          df <- data.frame(name = nm, stringsAsFactors = FALSE)
          df$value <- unname(list_for_df)
          df
        }
        
        geneSymbol3 <- list_to_df(geneSymbol)
        
        over2 <- cbind(over1,geneSymbol3)
        over3 <- over2[,-which(colnames(over2) %in% c("name"))]
        colnames(over3)[which(colnames(over3) %in% c("value"))] <- "geneSymbol"
        over3
      }
      
      over$bp <- convert2geneSymbol(over$bp)
      over$mf <- convert2geneSymbol(over$mf)
      over$cc <- convert2geneSymbol(over$cc)
      
      over_bp <- as_tibble(over$bp)
      over_mf <- as_tibble(over$mf)
      over_cc <- as_tibble(over$cc)
      
      over_bp <- over_bp[order(over_bp$BH.adjusted.p.value),]
      over_mf <- over_mf[order(over_mf$BH.adjusted.p.value),]
      over_cc <- over_cc[order(over_cc$BH.adjusted.p.value),]
      
      writeTibble <- function(tibble.input, output.file.name = tempfile())
      {
        if (!dir.exists(dirname(output.file.name)))
        {
          dir.create(dirname(output.file.name), recursive = TRUE)
        }
        flatten_list = function(x)
        {
          if (typeof(x) != "list")
          {
            return(x)
          }
          sapply(x, function(y) paste(y, collapse = " ; "))
        }
        tibble.input %>% mutate_all(funs(flatten_list)) %>% write.csv(output.file.name)
      }
      
      writeTibble(over_bp, output.file.name = file.path(output.file.dir,paste0(x_name,"_GO_BP.csv")))
      writeTibble(over_mf, output.file.name = file.path(output.file.dir,paste0(x_name,"_GO_MF.csv")))
      writeTibble(over_cc, output.file.name = file.path(output.file.dir,paste0(x_name,"_GO_CC.csv")))
    }
    
    #getGoAndPath(overlaps.anno.with.entrez.id,output.file.dir,x_name)
    getGoAndPath(overlaps.anno,output.file.dir,x_name)
    
  }else if(genome=="Hs"){
    
    library(EnsDb.Hsapiens.v75)
    #annoData<-toGRanges(EnsDb.Hsapiens.v75, feature="gene")
    
    dd.hs<-toGRanges(EnsDb.Hsapiens.v75)
    
    print(seqinfo(dd.hs))
    print(seqlevels(dd.hs))
    
    #print(seqlevels(dd.hs)[,1])
    
    #print(seqlevels(re.out[[1]])[,1])
    
    # seqlevels(dd.hs,force=TRUE) <- c("chr1","chr10","chr11","chr12","chr13",
    #                                           "chr14","chr15","chr16","chr17","chr18","chr19","chr2",
    #                                           "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY")
    
    #temp4=
    
    re.out.L<-lapply(1:length(re.out),function(u,re.out,dd.hs){
      
      x=re.out[[u]]
      x_name=names(re.out)[u]
      
      seqlevels(dd.hs,force=TRUE)<-seqinfo(x)@seqnames
      #print(seqinfo(re.out.trimmed))
      #print(seqlevels(re.out.trimmed))
      seqinfo(x)<-seqinfo(dd.hs)
      #GRCm38/mm10
      #dd<-toGRanges(EnsDb.Mmusculus.v79)
      #seqinfo(dd)
      #library(ensembldb)
      #library(GenomeInfoDb)
      seqlevelsStyle(x) <- seqlevelsStyle(dd.hs)
      re.out.trimmed<-trim(x, use.names=TRUE)
      overlaps.anno<-annoPeaks(re.out.trimmed,dd.hs)
      
      write.table(overlaps.anno,file=paste0(temp3,"/",x_name,"_annotation.txt"),row.names = FALSE,quote=FALSE,sep="\t")
    },re.out,dd.hs)
    
  }
  
},re.out,output.file.dir,genome)

}

# This function uses the overlapping peaks from findOverlapsOfPeak, and gets the 
# overlapping regions of two peak sets, and returns a GenomicRangesList including the
# following 3 GRanges object:
#  1: overlaping GRanges peaks
#  2: The peak1 GRanges peaks involoving overlapping  
#  3: The peak2 GRanges peaks involoving overlapping

getOverLap <- function(ol.ciLAD.XL.MEF_LAD,index_overlap,output.file.dir) {

  #index_overlap <- 1
  overlappedPeaks <- ol.ciLAD.XL.MEF_LAD$overlappingPeaks[[index_overlap]]
  print(names(overlappedPeaks))
  
  X_name <- names(ol.ciLAD.XL.MEF_LAD$overlappingPeaks)[index_overlap]
  
  XX_name <- gsub("-","_", X_name)
  
  XX_name <- gsub("///",",",XX_name)
  
  XXX_name <- unlist(strsplit(XX_name,","))
  
  print(XXX_name)
  
  colnames(overlappedPeaks)[which(colnames(overlappedPeaks)=="peaks1")] <- XXX_name[1]
  colnames(overlappedPeaks)[which(colnames(overlappedPeaks)=="peaks2")] <- XXX_name[2]
  
  overlap.peak.region <- apply(overlappedPeaks, 1, function(u){
    
    seqnames <- u["seqnames"]
    #start <- max(u["start"],u["start.1"])
    start <- max(u[names(u) %in% "start"])
    #end <- min(u["end"],u["end.1"])
    
    end <- min(u[names(u) %in% "end"])
    
    
    x <- data.frame(seqnames.olp=seqnames,start.olp=start,end.olp=end)      
    x
    
  })
  overlap.peak.region <- do.call(rbind.data.frame,overlap.peak.region)
  
  overlappedPeaks <- data.frame(overlappedPeaks,overlap.peak.region,check.names=F)
  
  overlappedPeaks <- data.frame(overlappedPeaks,overlappedBasePairs=as.numeric(as.character(overlappedPeaks$end.olp))
                                -as.numeric(as.character(overlappedPeaks$start.olp))+1,check.names=F)

  cat(colnames(overlappedPeaks),"\n")  
  
  output.file.name <- paste0(XXX_name[1],"-overlap-",XXX_name[2])
  
  output.file.name <- gsub("<","_", output.file.name)
  
  orderPeakAndOutPut(overlappedPeaks,output.file.dir,output.file.name,outHeader=T)
  
  #overlappedPeaks <- peaks.in.overlapping.region
  
  sample.index <- which(colnames(overlappedPeaks) == XXX_name[2])
  cat(sample.index,"\n")
  
  end <- sample.index-1
  peak1 <- unique(toGRanges(overlappedPeaks[,1:end],use.mcols=FALSE))
  
  overlap.index <- which(colnames(overlappedPeaks) == "seqnames.olp")
  end <- overlap.index-1
  peak2 <- unique(toGRanges(overlappedPeaks[,sample.index:end]))
  
  temp <- overlappedPeaks[,overlap.index:dim(overlappedPeaks)[2]]
  colnames(temp) <- c("seqnames","start","end","overlappedBasePairs")
  overlappingPeaks <- unique(toGRanges(temp))
  
  grl <- GenomicRangesList(peak1,peak2,overlappingPeaks)
  
  names(grl) <- c(XXX_name[1],XXX_name[2],X_name)
  grl
  
}

writeGrl2Bed <- function(grl,output.file.dir) {
  
  null <- lapply(1:length(grl), function(u,grl,output.file.dir){
    
     x <- names(grl)[u]
     y <- as.data.frame(grl[[u]])
     
     if(length(grep("///",x))==0){
       output.file.name  <- paste0(x,"_part.in.overlapping")
     }else{
       xx <- gsub("///","-overlaps-",x)
       output.file.name  <- paste0(xx,"_region")
     }
     
     orderPeakAndOutPut(y,output.file.dir,output.file.name,outHeader= FALSE)
     
  },grl,output.file.dir)
  
}

writeGrlInUnique2Bed <- function(grl,output.file.dir) {
  
  null <- lapply(1:length(grl), function(u,grl,output.file.dir){
    
    x <- names(grl)[u]
    y <- as.data.frame(grl[[u]])
    
    if(length(grep("///",x))==0){
      output.file.name  <- paste0(x,"_part.in.unique")
    }else{
      xx <- gsub("///","-overlaps-",x)
      output.file.name  <- paste0(xx,"_region")
    }
    
    orderPeakAndOutPut(y,output.file.dir,output.file.name,outHeader= FALSE)
    
  },grl,output.file.dir)
  
}

outGrl <- function(grl,output.file.dir) {

  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  null <- lapply(1:length(grl), function(u,grl,output.file.dir){
    x <- names(grl)[u]
    y <- as.data.frame(grl[[u]])
    output.file.name  <- x
    orderPeakAndOutPut(y,output.file.dir,output.file.name,outHeader= FALSE)
    
  },grl,output.file.dir)
  
}

# names(re.out)
# [1] "F121_9_all_1.5FC_F200_minp_adjPval-copy" "XL_MEF_1.5FC_minp_countF200_AvgSig-copy"

# ol.F121_9.XL_MEF <- findOverlapsOfPeaks(re.out[c(1,2)])

# ol.F121_9.XL_MEF.subsets <- getUniquePeaks(ol.F121_9.XL_MEF)

# output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_9_18_2018/F121_9-vs-XL_MEF"
# peaks.in.overlapping.region <- getOverLap(ol.F121_9.XL_MEF,1,output.file.dir)
# writeGrl2Bed(peaks.in.overlapping.region,output.file.dir)
# peaks.in.unique.regions.4.ol.F121_9.XL_MEF <- getPeaksInUnique(ol.F121_9.XL_MEF)
# writeGrlInUnique2Bed(peaks.in.unique.regions.4.ol.F121_9.XL_MEF,output.file.dir)

# names(nonXL.3.subsets) <- gsub("XL","nonXL",names(nonXL.3.subsets))
# 
# output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_9_18_2018/F121_9-vs-XL_MEF"
# null <- peaksToEnrichedGO(nonXL.3.subsets,output.file.dir,genome)

getPeaksInOverlap <- function(ol.ciLAD.XL.MEF_LAD,index_overlap) {
  
  overlappedPeaks <- ol.ciLAD.XL.MEF_LAD$overlappingPeaks[[index_overlap]]
  print(names(overlappedPeaks))
  
  X_name <- names(ol.ciLAD.XL.MEF_LAD$overlappingPeaks)[index_overlap]
  
  XX_name <- gsub("-","_", X_name)
  
  XX_name <- gsub("///",",",XX_name)
  
  XXX_name <- unlist(strsplit(XX_name,","))
  
  print(XXX_name)
  
  colnames(overlappedPeaks)[which(colnames(overlappedPeaks)=="peaks1")] <- XXX_name[1]
  colnames(overlappedPeaks)[which(colnames(overlappedPeaks)=="peaks2")] <- XXX_name[2]
  
  overlap.peak.region <- apply(overlappedPeaks, 1, function(u){
    
    seqnames <- u["seqnames"]
    start <- max(u[names(u) %in% "start"])
    end <- min(u[names(u) %in% "end"])
    
    x <- data.frame(seqnames.olp=seqnames,start.olp=start,end.olp=end)      
    x
    
  })
  overlap.peak.region <- do.call(rbind.data.frame,overlap.peak.region)
  
  overlappedPeaks <- data.frame(overlappedPeaks,overlap.peak.region,check.names=F)
  
  overlappedPeaks <- data.frame(overlappedPeaks,overlappedBasePairs=as.numeric(as.character(overlappedPeaks$end.olp))
                                -as.numeric(as.character(overlappedPeaks$start.olp))+1,check.names=F)
  

  sample.index <- which(colnames(overlappedPeaks) == XXX_name[2])

  end <- sample.index-1
  peak1 <- unique(toGRanges(overlappedPeaks[,1:end],use.mcols=FALSE))
  
  overlap.index <- which(colnames(overlappedPeaks) == "seqnames.olp")
  end <- overlap.index-1
  peak2 <- unique(toGRanges(overlappedPeaks[,sample.index:end]))
  
  temp <- overlappedPeaks[,overlap.index:dim(overlappedPeaks)[2]]
  colnames(temp) <- c("seqnames","start","end","overlappedBasePairs")
  overlappingPeaks <- unique(toGRanges(temp))
  
  grl <- GenomicRangesList(peak1,peak2,overlappingPeaks)
  
  names(grl) <- c(XXX_name[1],XXX_name[2],X_name)
  grl
  
}

getPeaksInAllOverlap <- function(ol.ciLAD.XL.MEF_LAD) {
  
  overlappedPeaks <- ol.ciLAD.XL.MEF_LAD$overlappingPeaks
  print(names(overlappedPeaks))
  
  s <- names(overlappedPeaks)
  
  zz <- unlist(lapply(1:length(s), function(u,ol.ciLAD.XL.MEF_LAD){
    
    z <- getPeaksInOverlap(ol.ciLAD.XL.MEF_LAD,u)
    z
  },ol.ciLAD.XL.MEF_LAD))
  
  zz
  
}

# XL.3.subsets.2 <- getUniquePeaks(ol.ciLAD.XL.MEF_LAD)
# ciLAD.3.subsets <- getUniquePeaks(ol.ciLAD.XL.MEF_LAD)
# LAD.3.subsets <- getUniquePeaks(ol.ciLAD.XL.MEF_LAD)

getUniquePeaks <- function(ol.ciLAD.XL.MEF_LAD) {
  
  ol.unique.0 <- getPeaksInUnique(ol.ciLAD.XL.MEF_LAD)
  
  zz <- getPeaksInAllOverlap(ol.ciLAD.XL.MEF_LAD)
    
  zz_name <- lapply(zz,function(u) {names(u)})
  
  print(zz_name)
  
  if(length(zz) > 1) {
    
  cat('Select the peak you want to overlap again:', '\n')
  m<-scan("",quiet=TRUE)
  
  cat("The peak names you selected are:",names(zz[[m[1]]])[m[2]],"\t",names(zz[[m[3]]])[m[4]],"\n")
  
  x <- paste0(toString((names(zz[[m[1]]])[m[2]])),"_in_",toString(names(zz[[m[1]]])[3]))
  y <- paste0(toString((names(zz[[m[3]]])[m[4]])),"_in_",toString(names(zz[[m[3]]])[3]))
  
  x <- gsub("///","_overlaps_",x)
  y <- gsub("///","_overlaps_",y)
  
  cat(x,"\t",y,"\n")
  
  w <- list(x=zz[[m[1]]][[m[2]]],y= zz[[m[3]]][[m[4]]])
  names(w) <- c(x,y)
 
  ol <- findOverlapsOfPeaks(w[c(1,2)],connectedPeaks = "keepAll")
  
  ol.unique.1 <- getPeaksInUnique(ol)
  
  XL.3.subsets <- c(ol.unique.0,ol.unique.1)
  } else{
        cat('Select the peak you want to include:', '\n')
        m<-scan("",quiet=TRUE)
        cat("The peak names you selected are:",names(zz[[m[1]]])[m[2]],"\n")
        x <- paste0(toString((names(zz[[m[1]]])[m[2]])),"_in_",toString(names(zz[[m[1]]])[3]))
        x <- gsub("///","_overlaps_",x)
        cat(x,"\n")
        w <- zz[[m[1]]][[m[2]]]
        XL.3.subsets <- c(ol.unique.0,w)
        names(XL.3.subsets)[2] <- x
        XL.3.subsets
  }
  
  XL.3.subsets
  
}

# fpkm.file <- "~/Aimin/DropboxUmass/NADfinder/BedFiles/GSM1621026_WT_fpkm.csv"
# fpkm.value <- dealWithRnaSeqFpkm(fpkm.file)

dealWithRnaSeqFpkm <- function(fpkm.file) {

    fpkm <- read.table(fpkm.file,sep = "\t",header = T)
  
    dat <- readLines(fpkm.file)
    dat <- strsplit(dat, ",")
  
  xx <- lapply(dat, function(u){
    
    n <- length(u)
    
    gene.name <- paste(u[1:n-1],collapse = ":")
    fpkm <- u[n]
    
    x <- data.frame(gene.name=gene.name,fpkm=fpkm)
    x 
       
  })
  xxx <- do.call(rbind,xx)
  xxx <- xxx[-1,]
  
#  class(fpkm$WT_FPKM)
  
#  which(is.na(fpkm$WT_FPKM))
  
  #fpkm.value <- as.numeric(as.character(xxx$fpkm))
  
  #hist(fpkm.value)
  #boxplot(fpkm.value)
  xxx
}

# output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_9_30_2018_reOrder_venn"  
# peak.index <- c(2,5)
# name <- c("H3K27me3","nonXL")
# getCount4Venn(re.out,peak.index,name,output.file.dir)

# peak.index <- c(3,5)
# name <- c("H3K9me3","nonXL")
# getCount4Venn(re.out,peak.index,name,output.file.dir)

# peak.index <- c(4,5)
# name <- c("PML","nonXL")
# getCount4Venn(re.out,peak.index,name,output.file.dir)

# peak.index <- c(2,7)
# name <- c("H3K27me3","XL")
# getCount4Venn(re.out,peak.index,name,output.file.dir)

# peak.index <- c(3,7)
# name <- c("H3K9me3","XL")
# getCount4Venn(re.out,peak.index,name,output.file.dir)

# peak.index <- c(4,7)
# name <- c("PML","XL")
# getCount4Venn(re.out,peak.index,name,output.file.dir)

# output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_10_1_2018_venn_with_ATAC"
# peak.index <- c(6,5)
# name <- c("nonXL","ATAC")
# getCount4Venn(re.out,peak.index,name,output.file.dir)

# peak.index <- c(8,5)
# name <- c("XL","ATAC")
# getCount4Venn(re.out,peak.index,name,output.file.dir)

# names(XL.nonXL.subset)
# names(re.out)

# NAD_class.ATAC <-  list(XL_uniq=XL.nonXL.subset[[1]],ATAC=re.out[[5]])
# peak.index <- c(1,2)
# name <- c("XL_uniq","ATAC")
# getCount4Venn(NAD_class.ATAC,peak.index,name,output.file.dir)

# NAD_class.ATAC <-  list(XL_ciLAD=XL.nonXL.subset[[2]],ATAC=re.out[[5]])
# peak.index <- c(1,2)
# name <- c("XL_ciLAD","ATAC")
# getCount4Venn(NAD_class.ATAC,peak.index,name,output.file.dir)

# NAD_class.ATAC <-  list(XL_LAD=XL.nonXL.subset[[3]],ATAC=re.out[[5]])
# peak.index <- c(1,2)
# name <- c("XL_LAD","ATAC")
# getCount4Venn(NAD_class.ATAC,peak.index,name,output.file.dir)

# NAD_class.ATAC <-  list(nonXL_uniq=XL.nonXL.subset[[4]],ATAC=re.out[[5]])
# peak.index <- c(1,2)
# name <- c("nonXL_uniq","ATAC")
# getCount4Venn(NAD_class.ATAC,peak.index,name,output.file.dir)

# NAD_class.ATAC <-  list(nonXL_ciLAD=XL.nonXL.subset[[5]],ATAC=re.out[[5]])
# peak.index <- c(1,2)
# name <- c("nonXL_ciLAD","ATAC")
# getCount4Venn(NAD_class.ATAC,peak.index,name,output.file.dir)

# NAD_class.ATAC <-  list(nonXL_LAD=XL.nonXL.subset[[6]],ATAC=re.out[[5]])
# peak.index <- c(1,2)
# name <- c("nonXL_LAD","ATAC")
# getCount4Venn(NAD_class.ATAC,peak.index,name,output.file.dir)




getCount4Venn <- function(re.out,peak.index,name,output.file.dir) {
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  grl <- GRangesList(re.out[peak.index])
  names(grl) <- name
  
  Z <- ChIPpeakAnno:::vennCounts(grl,length(peak.index),names=names(grl),by="base") 
  
  if(length(peak.index)==2){
    ZZ <- Z$venn_cnt[-which(row.names(Z$venn_cnt)=="00"),]
    labels = rev(name)
  }
  
  if(length(peak.index)==3){
    ZZ <- Z$venn_cnt[-which(row.names(Z$venn_cnt)=="000"),]
    labels = name
  }
  
  if(length(peak.index)==4){
    ZZ <- Z$venn_cnt[-which(row.names(Z$venn_cnt)=="0000"),]
    labels = name
  }
  
  y <- ZZ[,"Counts"]
  
  print(ZZ)
  
  names(y) <- row.names(ZZ)
  
  pdf(file = file.path(output.file.dir,paste0(paste(name,collapse = "-"),".pdf")))
  plot.new()
  colorfulVennPlot::plotVenn(y, labels, Colors=rainbow(18))
  dev.off()
  
}

# ol.ciLAD.XL.MEF_LAD <- findOverlapsOfPeaks(re.out[c(1,7,6)])
# XL.3.subsets <- getUniquePeaks(ol.ciLAD.XL.MEF_LAD)
# XL.3.subsets.anno <- getAnnotatedGene(XL.3.subsets,"Mm")

# ol.ciLAD.non_XL.MEF_LAD <- findOverlapsOfPeaks(re.out[c(1,5,6)])
# nonXL.3.subsets <- getUniquePeaks(ol.ciLAD.non_XL.MEF_LAD)
# nonXL.3.subsets.anno <- getAnnotatedGene(nonXL.3.subsets,"Mm")

getAnnotatedGene <- function(re.out,genome){
  
  y <- lapply(1:length(re.out),function(u,re.out,genome){
    
    x=re.out[[u]]
    x_name=names(re.out)[u]
    
    if(genome=="Mm"){
      
      dd.GRCm39.mm10<-toGRanges(EnsDb.Mmusculus.v75)
      
      overlaps.trimmed<-trim(x,use.names=TRUE)
      overlaps.anno<-annoPeaks(overlaps.trimmed,dd.GRCm39.mm10)
     
    }
    overlaps.anno
  },re.out,genome)

  names(y) <- names(re.out)
  y
  
}

# YYY1 <- getFPKM4DiffSet(XL.3.subsets.anno,fpkm.value)
# YYY2 <- getFPKM4DiffSet(nonXL.3.subsets.anno,fpkm.value)
# YYY <- rbind(YYY1,YYY2)

getFPKM4DiffSet <- function(XL.3.subsets.anno,fpkm.value) {
  YY <- lapply(1:length(XL.3.subsets.anno), function(u,XL.3.subsets.anno,fpkm.value){
  
    index.matched.2.fpkm <- match(XL.3.subsets.anno[[u]]$gene_name,fpkm.value$gene.name)
    index1 <- index.matched.2.fpkm[which(!is.na(index.matched.2.fpkm))]
    fpkm.value.4.oneset <- fpkm.value[index1,]
    set.name <- names(XL.3.subsets.anno)[u]
    y <- cbind.data.frame(rep(set.name,dim(fpkm.value.4.oneset)[1]),fpkm.value.4.oneset)
    colnames(y) <- c("SetName","GeneName","FPKM")
    y
  },XL.3.subsets.anno,fpkm.value)
  
  YYY <- do.call(rbind,YY)
  YYY
}

# getBoxPlot4FPKMOfSubsetPeaks(YYY,output.file.dir)
# output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_10_2_2018_boxplot_include_whole_genome"
# YYYY <- getBoxPlot4FPKMOfSubsetPeaks(YYY,fpkm.value,output.file.dir)
 
getBoxPlot4FPKMOfSubsetPeaks <- function(YYY,fpkm.value,output.file.dir){
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
YYY.set.name <- unique(as.character(YYY$SetName))

index <- which(levels(YYY$SetName)==YYY.set.name[1])
levels(YYY$SetName)[index] <- "XL_Uniq"

index <- which(levels(YYY$SetName)==YYY.set.name[2])
levels(YYY$SetName)[index] <- "XL_ciLAD"

index <- which(levels(YYY$SetName)==YYY.set.name[3])
levels(YYY$SetName)[index] <- "XL_LAD"

index <- which(levels(YYY$SetName)==YYY.set.name[4])
levels(YYY$SetName)[index] <- "nonXL_Uniq"

index <- which(levels(YYY$SetName)==YYY.set.name[5])
levels(YYY$SetName)[index] <- "nonXL_ciLAD"

index <- which(levels(YYY$SetName)==YYY.set.name[6])
levels(YYY$SetName)[index] <- "nonXL_LAD"

fpkm.value.ref <- data.frame(SetName=rep("wholeGenome",dim(fpkm.value)[1]),GeneName=fpkm.value$gene.name,FPKM=fpkm.value$fpkm)

NAD.fpkm.value <- fpkm.value.ref[which(fpkm.value.ref$GeneName %in% YYY$GeneName),]
NAD.fpkm.value$SetName <- "NAD"
nonNAD.fpkm.value <- fpkm.value.ref[-which(fpkm.value.ref$GeneName %in% YYY$GeneName),]
nonNAD.fpkm.value$SetName <- "nonNAD"

YYYY <- rbind(YYY,NAD.fpkm.value,nonNAD.fpkm.value,fpkm.value.ref)

png(file = file.path(output.file.dir,paste0(paste(levels(YYYY$SetName),collapse = "-"),".png")),width = 1500, height = 480)
boxplot(log10(as.numeric(as.character(YYYY$FPKM))+1)~YYYY$SetName,ylab = "log10(FPKM+1)")
dev.off()

YYYY

}

# ProcessUsingLOLA

# ctcf.input.file <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/CtcfCollection"
# selected.peak.index <- c(4,7:10,13)
# dbPath = "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/CtcfCollection/mm10"
# regionDB <- generateRegionDB(ctcf.input.file,selected.peak.index,dpPath)

generateRegionDB <- function(ctcf.input.file,selected.peak.index,dpPath) {
  
  ctcf.input.file.name <- list.files(ctcf.input.file,all.files = TRUE,full.names = TRUE,recursive = FALSE,include.dirs = TRUE)
  ctcf.input.file.name.1 <- ctcf.input.file.name[selected.peak.index]
  ctcf.bed.collection.dir <- file.path(dbPath,"/ucsc_example/regions")

  if(!dir.exists(ctcf.bed.collection.dir)){dir.create(ctcf.bed.collection.dir,recursive = TRUE)}
  
  null <- lapply(1:length(ctcf.input.file.name.1), function(u,ctcf.input.file.name.1,ctcf.bed.collection.dir){
    
    x <- tools::file_path_sans_ext(basename(ctcf.input.file.name.1[u]))
    
    cmd = paste("grep -v track",ctcf.input.file.name.1[u],"| awk -v OFS='\t' -F'\t' '{print $1,$2,$3}'",">",file.path(ctcf.bed.collection.dir,paste0(x,".bed")),collapse= " ")
    system(cmd)  
    
  },ctcf.input.file.name.1,ctcf.bed.collection.dir)
  
  #dbPath = "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/CtcfCollection/mm10"
  
  regionDB = loadRegionDB(dbLocation=dbPath)
  regionDB
}

# output.file.dir <-  "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_9_28_2018"
# dbPath = "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/CtcfCollection/mm10"
# regionDB = loadRegionDB(dbLocation=dbPath)
# XL.nonXL.subset <- c(XL.3.subsets,nonXL.3.subsets)
# getOverlapWithOther(XL.nonXL.subset,regionDB,output.file.dir)

# select.query.peak.index <- c(1,2,3)
# output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_10_1_2018_WithCTCF"


# output.file.dir <-  "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_10_3_2018_include_Non_CTCF"
# XL.nonXL.subset <- c(XL.3.subsets,nonXL.3.subsets)
# 



getOverlapWithOther <- function(XL.nonXL.subset,select.query.peak.index= NULL,regionDB,output.file.dir) {
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  grl <- lapply(1:length(XL.nonXL.subset), function(u,XL.nonXL.subset){
     
    gr <- XL.nonXL.subset[[u]]
    mcols(gr) <- NULL
    gr
      
  },XL.nonXL.subset)
  names(grl) <- names(XL.nonXL.subset)
  
  grl2 <- GRangesList(grl)
  
  if(is.null(select.query.peak.index)){
    
    grl3 <- c(regionDB$regionGRL,grl2)
    
    userUniverse <- buildRestrictedUniverse(grl3)
    
  locResults = runLOLA(grl2,userUniverse,regionDB,cores=4)
  
  userSetName <- unique(locResults$userSet)
  
  locResults[which(locResults$userSet==userSetName[1]),]$userSet <- "XL_LAD"
  locResults[which(locResults$userSet==userSetName[2]),]$userSet <- "nonXL_LAD"
  locResults[which(locResults$userSet==userSetName[3]),]$userSet <- "XL_ciLAD"
  locResults[which(locResults$userSet==userSetName[4]),]$userSet <- "nonXL_ciLAD"
  locResults[which(locResults$userSet==userSetName[5]),]$userSet <- "XL"
  locResults[which(locResults$userSet==userSetName[6]),]$userSet <- "nonXL"
  
  write.table(locResults,
              file=file.path(output.file.dir,"OverlapWith.other.CTCF.data.set.txt"),sep="\t",
              quote = FALSE,row.names = FALSE,col.names = TRUE)
  }else{
    
    select.query.peak.index <- c(1,2,3)
    
    grl3 <- c(regionDB$regionGRL,grl2[select.query.peak.index])
    
    userUniverse <- buildRestrictedUniverse(grl3)
    
    
    locResults = runLOLA(grl2[select.query.peak.index],userUniverse,regionDB,cores=4)
    
    userSetName <- unique(locResults$userSet)
    
    locResults[which(locResults$userSet==userSetName[1]),]$userSet <- "XL_LAD"
    locResults[which(locResults$userSet==userSetName[2]),]$userSet <- "XL_ciLAD"
    locResults[which(locResults$userSet==userSetName[3]),]$userSet <- "XL"

    write.table(locResults,
                file=file.path(output.file.dir,"XL_OverlapWith.other.CTCF.data.set.txt"),sep="\t",
                quote = FALSE,row.names = FALSE,col.names = TRUE)
    
    select.query.peak.index <- c(4,5,6)
    
    grl3 <- c(regionDB$regionGRL,grl2[select.query.peak.index])
    
    userUniverse <- buildRestrictedUniverse(grl3)
    
    locResults = runLOLA(grl2[select.query.peak.index],userUniverse,regionDB,cores=4)
    
    userSetName <- unique(locResults$userSet)
    
    locResults[which(locResults$userSet==userSetName[1]),]$userSet <- "nonXL_LAD"
    locResults[which(locResults$userSet==userSetName[2]),]$userSet <- "nonXL_ciLAD"
    locResults[which(locResults$userSet==userSetName[3]),]$userSet <- "nonXL"
    
    write.table(locResults,
                file=file.path(output.file.dir,"nonXL_OverlapWith.other.CTCF.data.set.txt"),sep="\t",
                quote = FALSE,row.names = FALSE,col.names = TRUE)
    
  }

  locResults
}

# Make venns for the following:
# output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_10_1_2018_venn_with_modifications"

# nonXL, H3K27me3, ciLAD
# peak.index <- c(5,2,1)
# name <- c("nonXL","H3K27me3","ciLAD")
# getCount4Venn(re.out,peak.index,name,output.file.dir)

# nonXL, H3K27me3, LAD
# peak.index <- c(5,2,6)
# name <- c("nonXL","H3K27me3","LAD")
# getCount4Venn(re.out,peak.index,name,output.file.dir)

# nonXL, H3K9me3, ciLAD
# peak.index <- c(5,3,1)
# name <- c("nonXL","H3K9me3","ciLAD")
# getCount4Venn(re.out,peak.index,name,output.file.dir)

# nonXL, H3K9me3, LAD
# peak.index <- c(5,3,6)
# name <- c("nonXL","H3K9me3", "LAD")
# getCount4Venn(re.out,peak.index,name,output.file.dir)


# input.bed.dir <- "/Users/aiminyan/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/bedFilesforDensityPlots/byClass"

# input.bw.path <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/H3K27me3_bigwig"

# output.file.dir <- "~/Aimin/umw_nathan_lawson/Aimin/H3K27me3_cvg_3"
# cvglists.l <- overLapWithOtherFeatures(input.file.dir,input.bw.path,output.file.dir)
# 
# output.file.dir <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/Aimin/DensityPlot"

# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/test"

# input.bed.dir <- "/Users/aiminyan/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/bedFilesforDensityPlots/byClass"
# input.bw.path <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/H3K27me3_bigwig"
# dd <- 3000
# dd <- 10000
# tfType <- "H3K27me3"
# class_pattern <- "tss"
# tssType <- c("Class1","Class2","Class3","Endo")
# output.file.dir <- paste0("~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/H3K27me3_",class_pattern)
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType)

# dd <- 3000
# dd <- 10000
# class_pattern <- "enh"
# output.file.dir <- paste0("~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/H3K27me3_",class_pattern)
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern)

# dd <- 3000
# dd <- 10000
# class_pattern <- "DHS"
# output.file.dir <- paste0("~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/H3K27me3_",class_pattern)
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern)


# ArteryTSS.bed
# VeinTSS.bed
# CommonTSS.bed
# ArteryEnh.bed
# VeinEnh.bed
# CommonEnh.bed

# input.bed.dir <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/bedFilesforDensityPlots/AllElementsByCellType/"
# input.bw.path <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/H3K27me3_bigwig"
# dd <- 3000
# dd <- 10000
# class_pattern <- "TSS"
# tssType <- c("Artery","Common","Vein")
# output.file.dir <- paste0("~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/H3K27me3_",paste0("Big_",class_pattern))
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType)

# dd <- 3000
# dd <- 10000
# class_pattern <- "Enh"
# tssType <- c("Artery","Common","Vein")
# output.file.dir <- paste0("~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/H3K27me3_",paste0("Big_",class_pattern))
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType)


# input.bed.dir <- "/Users/aiminyan/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/bedFilesforDensityPlots/byClass"
# input.bw.path <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/NR2F2kd_bigwig"
# dd <- 3000
# dd <- 10000
# tfType <- "TRC-NS"
# class_pattern <- "tss"
# tssType <- c("Class1","Class2","Class3","Endo")
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/NR2F2"
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

# dd <- 3000
# dd <- 10000
# tfType <- "TRC-NS"
# class_pattern <- "enh"
# tssType <- c("Class1","Class2","Class3","Endo")
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/NR2F2"
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

# dd <- 3000
# dd <- 10000
# tfType <- "TRC-NS"
# class_pattern <- "DHS"
# tssType <- c("Class1","Class2","Class3","Endo")
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/NR2F2"
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

# dd <- 3000
# dd <- 10000
# tfType <- "shNR2F2"
# class_pattern <- "enh"
# tssType <- c("Class1","Class2","Class3","Endo")
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/NR2F2"
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

# dd <- 3000
# dd <- 10000
# tfType <- "shNR2F2"
# class_pattern <- "DHS"
# tssType <- c("Class1","Class2","Class3","Endo")
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/NR2F2"
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

# tfType <- "shNR2F2"
# class_pattern <- "tss"
# tssType <- c("Class1","Class2","Class3","Endo")
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/NR2F2"
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

# input.bed.dir <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/bedFilesforDensityPlots/AllElementsByCellType/"
# input.bw.path <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/NR2F2kd_bigwig"
# dd <- 3000
# dd <- 10000
# tfType <- "TRC-NS"

# class_pattern <- "TSS"
# tssType <- c("Artery","Common","Vein")
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/NR2F2"
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

# class_pattern <- "Enh"
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/NR2F2"
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

# dd <- 3000
# dd <- 10000
# tfType <- "shNR2F2"

# class_pattern <- "TSS"
# tssType <- c("Artery","Common","Vein")
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/NR2F2"
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

# class_pattern <- "Enh"
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/NR2F2"
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

# input.bed.dir <- "/Users/aiminyan/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/bedFilesforDensityPlots/byClass"
# input.bw.path <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/Aimin/EZH2bam_ENCODE_bw"
# 
# dd <- 3000
# dd <- 10000
# tfType <- "EZH2"
# class_pattern <- "tss"
# class_pattern <- "enh"
# class_pattern <- "DHS"
# tssType <- c("Class1","Class2","Class3","Endo")


# input.bed.dir <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/bedFilesforDensityPlots/AllElementsByCellType/"
#
# class_pattern <- "TSS"
# tssType <- c("Artery","Common","Vein")

#
# class_pattern <- "Enh"
# tssType <- c("Artery","Common","Vein")

# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/EZH2"
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

# input.bed.dir <- "/Users/aiminyan/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/bedFilesforDensityPlots/byClass"
# input.bw.path <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/Aimin/HDAC1cr_bw"
# 
# output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/HDAC1"
# dd <- 3000
# dd <- 10000
# tfType <- "HDAC1"
# class_pattern <- "tss"
# class_pattern <- "enh"
# class_pattern <- "DHS"
# tssType <- c("Class1","Class2","Class3","Endo")
# null <- overLapWithOtherFeatures(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType)

overLapWithOtherFeatures <- function(input.bed.dir,input.bw.path,output.file.dir,dd,class_pattern,tssType,tfType) {
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  file.name.4 <- list.files(input.bed.dir,pattern=".bed$",all.files = TRUE,full.names = TRUE,recursive = FALSE,include.dirs = TRUE)
  
  data <- sapply(file.name.4, toGRanges, format="BED")
  
  names(data) <- gsub(".bed", "", basename(names(data)))
  
  #tss <- c("class1tss_2KplMin","class2tss_2KplMin","class3tss_2KplMin","EndoRefTSS_2KplMin")
  
  #tss <- names(data)[grep("tss",names(data),ignore.case = T)]
  
  tss <- names(data)[grep(class_pattern,names(data),ignore.case = T)]
  
  tss.bed <- data[which(names(data) %in% tss)]

  feature.center <- lapply(1:length(tss.bed), function(u,tss.bed,dd){
    x <- tss.bed[[u]]
    x
  },tss.bed,dd)
  
  names(feature.center) <- names(tss.bed)
  
  # get the regions of around feature middle +/- dd 
  feature.recentered.l <- lapply(1:length(tss.bed), function(u,tss.bed,dd){
    
    features <- tss.bed[[u]]
    
    wid <- width(features)
    feature.recentered <- feature.center <- features
    
    start(feature.center) <- start(features) + floor(wid/2)
    
    width(feature.center) <- 1
    
    start(feature.recentered) <- start(feature.center) - dd
    end(feature.recentered) <- end(feature.center) + dd
    
    feature.recentered
    
  },tss.bed,dd)
  
  names(feature.recentered.l) <- names(tss.bed)
  
  # feature.test <- tss.bed[[1]]
  # wid.test <- width(feature.test)
  # feature.recentered.test <- feature.center.test <- feature.test
  # start(feature.center.test) <- start(feature.test) + floor(wid.test/2)
  # width(feature.center.test) <- 1
  # start(feature.recentered.test) <- start(feature.center.test) - dd
  # end(feature.recentered.test) <- end(feature.center.test) + dd
  
  cvglists.l <- lapply(1:length(feature.recentered.l), function(u,feature.recentered.l,input.bw.path){
    
    files <- dir(input.bw.path, "bw")
    if(.Platform$OS.type != "windows"){
      cvglists <- sapply(file.path(input.bw.path,files), import, 
                         format="BigWig", 
                         which=feature.recentered.l[[u]], 
                         as="RleList")
    }else{## rtracklayer can not import bigWig files on Windows
      load(file.path(input.bw.path, "cvglist.rds"))
    }
    
  },feature.recentered.l,input.bw.path)
  
  
  grep(tfType,names(unlist(cvglists.l)))
  
  
  #x_name_0 <- tools::file_path_sans_ext(basename(names(cvglists.l[[1]])))  
  #x_name_start <- str_locate_all(x_name_0,"_")[[1]][2,2] + 1
  #x_name_end <- str_locate_all(x_name_0,"_")[[1]][3,2] - 1
  #x_name <- paste0(unique(str_sub(x_name_0,x_name_start,x_name_end)),"_",dd,"_",paste(tss,collapse = "_"))
  
  x_name <- paste0(tfType,"_",dd,"_",paste(tss,collapse = "_"))
  
  saveRDS(cvglists.l,file = file.path(output.file.dir,paste0(x_name,"_cvglist.rds")))
  
  
  
  #lapply(1:length(tss.bed),function(u){
  
  #Class1 <- (cvglists.ul[[1]]+cvglists.ul[[2]])/2
  
  #Class2 <- (cvglists.ul[[3]]+cvglists.ul[[4]])/2
  
  #Class3 <- (cvglists.ul[[5]]+cvglists.ul[[6]])/2
  
  #Endo <- (cvglists.ul[[7]]+cvglists.ul[[8]])/2
  
  #}
  
  #cvglists.ave <- list(Class1=Class1,Class2=Class2,Class3=Class3, Endo=Endo)
  
  cvglists.l <- lapply(1:length(cvglists.l), function(u,cvglists.l){
    
    z <- cvglists.l[[u]][grep(tfType,names(cvglists.l[[u]]))]
    z
  },cvglists.l)
  
  #names(cvglists.l[[4]])
  
  cvglists.ul <- unlist(cvglists.l)
  
  #cvglists.ul <- cvglists.ul[grep(tfType,names(cvglists.ul))]
  
  #grep(tfType,names(unlist(cvglists.l)))
  
  n <- length(cvglists.ul)/length(tss.bed)
  
  if(n==2){
  cvglists.ave <- lapply(1:length(tss.bed),function(u,cvglists.ul){
    
    Class1 <- (cvglists.ul[[2*u-1]]+cvglists.ul[[2*u]])/n
    
    Class1
    
  },cvglists.ul)
  }
  
  if(n==3){
    cvglists.ave <- lapply(1:length(tss.bed),function(u,cvglists.ul){
      
      Class1 <- (cvglists.ul[[3*u-2]]+cvglists.ul[[3*u-1]]+cvglists.ul[[3*u]])/n
      
      Class1
      
    },cvglists.ul)
  }

  #cvglists.ave <- list(Class1=Class1,Class2=Class2,Class3=Class3, Endo=Endo)
  
  names(cvglists.ave) <- names(tss.bed)
  
  
  
  
  # sig1 <- featureAlignedSignal(cvglists.ave, feature.center[[1]], upstream=dd, downstream=dd) 
  # sigs.log2 <- lapply(sig1, function(.ele) log2(.ele+1))   
  # out.1 <- featureAlignedDistribution(sigs.log2, 
  #                                   feature.center[[1]],upstream=dd, downstream=dd,
  #                                   zeroAt=.5,type="l", 
  #                                   ylab="Averaged coverage")
  # 
  # 
  # sig1 <- featureAlignedSignal(cvglists.ave, feature.center[[2]], upstream=dd, downstream=dd) 
  # sigs.log2 <- lapply(sig1, function(.ele) log2(.ele+1))
  # out.2 <- featureAlignedDistribution(sigs.log2, 
  #                                   feature.center[[2]],upstream=dd, downstream=dd,
  #                                   zeroAt=.5,type="l", 
  #                                   ylab="Averaged coverage",add=TRUE)
  # 
  # sig1 <- featureAlignedSignal(cvglists.ave, feature.center[[3]], upstream=dd, downstream=dd) 
  # sigs.log2 <- lapply(sig1, function(.ele) log2(.ele+1))
  # out.3 <- featureAlignedDistribution(sigs.log2, 
  #                                   feature.center[[3]],upstream=dd, downstream=dd,
  #                                   zeroAt=.5,type="l", 
  #                                   ylab="Averaged coverage",add=TRUE)
  # 
  # 
  # sig1 <- featureAlignedSignal(cvglists.ave, feature.center[[4]], upstream=dd, downstream=dd) 
  # sigs.log2 <- lapply(sig1, function(.ele) log2(.ele+1))
  # out.4 <- featureAlignedDistribution(sigs.log2, 
  #                                   feature.center[[4]],upstream=dd, downstream=dd,
  #                                   zeroAt=.5,type="l", 
  #                                   ylab="Averaged coverage",add=TRUE)
  # 
  # out <- cbind(out.1$density[,1],out.2$density[,2],out.3$density[,3],out.4$density[,4])
  
  XX <- lapply(1:length(feature.center),function(u,feature.center){
    
    sig1 <- featureAlignedSignal(cvglists.ave, feature.center[[u]], upstream=dd, downstream=dd) 
    sigs.log2 <- lapply(sig1, function(.ele) log2(.ele+1))
    out.4 <- featureAlignedDistribution(sigs.log2, 
                                       feature.center[[u]],upstream=dd, downstream=dd,
                                       zeroAt=.5,type="l", 
                                       ylab="Averaged coverage",add=TRUE)
    
    x <- list(density=out.4$density[,u],grWidAt=out.4$grWidAt,grWidLab=out.4$grWidLab)
    x 
  },feature.center)
  
  den.l <-lapply(1:length(XX),function(u,XX){
    Z <- XX[[u]]$density
    Z
  },XX)
  out <- do.call(cbind,den.l)
  
  grWidAt.l <-lapply(1:length(XX),function(u,XX){
    Z <- XX[[u]]$grWidAt
    Z
  },XX)
  grWidAt <- do.call(cbind,grWidAt.l)
  
  grWidLab.l <-lapply(1:length(XX),function(u,XX){
    Z <- XX[[u]]$grWidLab
    Z
  },XX)
  grWidLab <- do.call(cbind,grWidLab.l)
  
  
  head(out)
  colnames(out) <- tssType

  drowDensity <- function(density,grWidAt,grWidLab,...) {
    dots <- list(...)
    matplot(density, ..., xaxt="n")
    axis(1, at = grWidAt, labels = grWidLab)
    lty <- if(!is.null(dots$lty)) dots$lty else 1:5
    lwd <- if(!is.null(dots$lwd)) dots$lwd else 1
    col <- if(!is.null(dots$col)) dots$col else 1:6
    legend("topright", legend=colnames(density), col=col,
           lty=lty, lwd=lwd)
  }

  jpeg(file.path(output.file.dir,paste0(x_name,"_density.jpeg")))
  drowDensity(out,grWidAt[,1],grWidLab[,1],type="l", lty=1,
              ylab="log2(averaged coverage)",xlab=paste0("Distance to ",class_pattern,"(bp)"),main=tfType)
  dev.off()
  
  x.pos <- seq(-dd,dd,length.out = dim(out)[1])
  out.m <- data.frame(x.pos=x.pos,out)
  
  outLong <- melt(data          = out.m,
                  id.vars       = c("x.pos"),
                  measure.vars  = tssType,
                  variable.name = "TSSByClass",
                  value.name    = "value")
  
  ggplot(outLong, aes(x.pos, value,colour=TSSByClass)) + 
    geom_line() + labs(x="Distance to TSS(bp)",y="log2(averaged coverage)",title=tfType) +
    theme(plot.title=element_text(hjust=0.5)) + scale_x_continuous(name=paste0("Distance to ",class_pattern,"(bp)"), limits=c(-dd, dd))
 
  write.table(out.m,file = file.path(output.file.dir,paste0(x_name,"_Density_matrix_short",".txt")),
              append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T)
  
  write.table(outLong,file = file.path(output.file.dir,paste0(x_name,"_Density_matrix_long",".txt")),
              append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T)
  
  null <- lapply(1:length(cvglists.l), function(u,cvglists.l,feature.center,dd,n,tfType){
    #u <- 1
  sig1.4.heatmap <- featureAlignedSignal(cvglists.l[[u]], feature.center[[u]], upstream=dd, downstream=dd) 
  sig1.4.heatmap.log2 <- lapply(sig1.4.heatmap, function(.ele) log2(.ele+1))
  names(sig1.4.heatmap.log2) <- paste0(tfType,"_",seq(1,n))
  
  #u <- 1
  #x_name_0 <- tools::file_path_sans_ext(basename(names(cvglists.l[[u]])))  
  #x_name_start <- str_locate_all(x_name_0,"_")[[1]][2,2] + 1
  #x_name_end <- str_locate_all(x_name_0,"_")[[1]][3,2] - 1

  x_name <- paste0(tfType,"_",names(feature.center)[u])
  
  jpeg(file.path(output.file.dir,paste0(x_name,"_",dd,"_heatmap.jpeg")))
  featureAlignedHeatmap(sig1.4.heatmap.log2, feature.center[[u]],upstream=dd, downstream=dd,zeroAt=.5,res=300)
  dev.off()
  
  },cvglists.l,feature.center,dd,n,tfType)
  
}

getNonCtcfBed <- function(re.out, BSgenome.Mmusculus.UCSC.mm10) {
  #paul.dir <- "/project/umw_paul_kaufman"
  pau.dir <- "~/Aimin/ProjectAtCluster/umw_paul_kaufman/Aimin"
  
  if(!dir.exists(pau.dir)){dir.create(pau.dir,recursive = TRUE)}
  
  null <- lapply(1:length(re.out), function(u,re.out){
    export(re.out[[u]],file.path(pau.dir,paste0(names(re.out)[u],".bed")),format="BED")},re.out)
  
  gr.whole.genome <- GRanges(seqinfo(BSgenome.Mmusculus.UCSC.mm10))
  
  bed.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/CtcfCollection"
  export(gr.whole.genome,file.path(bed.dir,"WholeGenome.bed"),format="BED")
  
  ref.bed <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/CtcfCollection/WholeGenome.bed"
  exclude.bed.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/CtcfCollection/mm10/ucsc_example/regions"
  
  exclude.bed.file.name <- list.files(exclude.bed.dir,all.files = TRUE,full.names = TRUE,recursive = FALSE,include.dirs = TRUE)
  
  exclude.bed <- paste(exclude.bed.file.name[4:7],collapse = " ")
  
  output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/CtcfCollection"
  cmd = paste("bedops -d",ref.bed,exclude.bed,">",file.path(output.file.dir,"NonCTCF.bed"))
  
  system(cmd)
}

# bed.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Results_10_4_2018_gene_density"

# bed.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Results_10_7_2018_gene_density"
# getGeneDensityPlot(bed.dir)

getGeneDensityPlot <- function(bed.dir) {
  
  if(!dir.exists(bed.dir)){dir.create(bed.dir,recursive = TRUE)}
  # whole genome in MB
  whole.genome.MB <- sum(width(gr.whole.genome))/1000000
  
  names(XL.nonXL.subset) <- c("XL_Uniq","XL_ciLAD","XL_LAD","nonXL_Uniq","nonXL_ciLAD","nonXL_LAD")
  
  null <- lapply(1:length(XL.nonXL.subset), function(u,XL.nonXL.subset){
    
    x <- names(XL.nonXL.subset)[u] 
    export(XL.nonXL.subset[[u]],file.path(bed.dir,paste0(x,".bed")),format="BED")
  
    },XL.nonXL.subset)
  
  xl.nad.bed <- file.path(bed.dir,paste0(c("XL_Uniq","XL_ciLAD","XL_LAD"),".bed"))
  cmd = paste("bedops -u",paste(xl.nad.bed,collapse = " "),">",file.path(bed.dir,"XLNAD.bed"))
  system(cmd)
  cmd = paste("bedtools sort -i",file.path(bed.dir,"XLNAD.bed"),">",file.path(bed.dir,"XLNAD_sorted.bed"))
  system(cmd) 
  cmd = paste("bedtools merge -i",file.path(bed.dir,"XLNAD_sorted.bed"),">",file.path(bed.dir,"XLNAD_merged.bed"))
  system(cmd) 
  cmd = paste("bedops -d",ref.bed,file.path(bed.dir,"XLNAD_merged.bed"),">",file.path(bed.dir,"XL_nonNAD.bed"))
  system(cmd)
  
  
  peaks=read.table(file.path(bed.dir,"XLNAD_merged.bed"))
  colnames(peaks)[1:3]= c("chr","start","end")
  XLNAD=toGRanges(peaks)
  
  peaks=read.table(file.path(bed.dir,"XL_nonNAD.bed"))
  colnames(peaks)[1:3]= c("chr","start","end")
  XL_nonNAD=toGRanges(peaks)
  
  
  
  nonxl.nad.bed <- file.path(bed.dir,paste0(c("nonXL_Uniq","nonXL_ciLAD","nonXL_LAD"),".bed"))
  cmd = paste("bedops -u",paste(nonxl.nad.bed,collapse = " "),">",file.path(bed.dir,"nonXLNAD.bed"))
  system(cmd)
  cmd = paste("bedtools sort -i",file.path(bed.dir,"nonXLNAD.bed"),">",file.path(bed.dir,"nonXLNAD_sorted.bed"))
  system(cmd) 
  cmd = paste("bedtools merge -i",file.path(bed.dir,"nonXLNAD_sorted.bed"),">",file.path(bed.dir,"nonXLNAD_merged.bed"))
  system(cmd) 
  cmd = paste("bedops -d",ref.bed,file.path(bed.dir,"nonXLNAD_merged.bed"),">",file.path(bed.dir,"nonXL_nonNAD.bed"))
  system(cmd)
  
  peaks=read.table(file.path(bed.dir,"nonXLNAD_merged.bed"))
  colnames(peaks)[1:3]= c("chr","start","end")
  nonXLNAD=toGRanges(peaks)
  
  peaks=read.table(file.path(bed.dir,"nonXL_nonNAD.bed"))
  colnames(peaks)[1:3]= c("chr","start","end")
  nonXL_nonNAD=toGRanges(peaks)
  
  XLNAD.MB <- sum(width(XLNAD))/1000000
  XL_nonNAD.MB <- sum(width(XL_nonNAD))/1000000
  
  nonXLNAD.MB <- sum(width(nonXLNAD))/1000000
  nonXL_nonNAD.MB <- sum(width(nonXL_nonNAD))/1000000
  
  
  # export(re.out[[6]],file.path(bed.dir,"nonXL.bed"),format="BED")
  # export(re.out[[8]],file.path(bed.dir,"XL.bed"),format="BED")
  # 
  # nad.bed <- paste(file.path(bed.dir,"nonXL.bed"),file.path(bed.dir,"XL.bed"),sep= " ")
  # cmd = paste("bedops -u",nad.bed,">",file.path(bed.dir,"NAD.bed"))
  # system(cmd)
  # cmd = paste("bedtools sort -i",file.path(bed.dir,"NAD.bed"),">",file.path(bed.dir,"NAD_sorted.bed"))
  # system(cmd) 
  # cmd = paste("bedtools merge -i",file.path(bed.dir,"NAD_sorted.bed"),">",file.path(bed.dir,"NAD_merged.bed"))
  # system(cmd) 
  # peaks=read.table(file.path(bed.dir,"NAD_merged.bed"))
  # colnames(peaks)[1:3]= c("chr","start","end")
  # NAD=toGRanges(peaks)
  # 
  # # NAD.before <- unique(c(re.out[[6]],re.out[[8]]))
  # 
  # XLNAD <- re.out[[8]]
  # nonXLNAD <- re.out[[6]]
  #   
  # # export(NAD,file.path(bed.dir,"NAD.bed"),format="BED")
  # 
  # 
  # cmd = paste("bedops -d",ref.bed,file.path(bed.dir,"NAD_merged.bed"),">",file.path(bed.dir,"nonNAD.bed"))
  # 
  # system(cmd)
  # 
  # peaks=read.table(file.path(bed.dir,"nonNAD.bed"))
  # 
  # colnames(peaks)[1:3]= c("chr","start","end")
  # nonNAD=toGRanges(peaks)
  # 
  # #mcols(NAD) <- NULL
  # #nonNAD <- setdiff(gr.whole.genome,NAD)
  # 
  # XLNAD.MB <- sum(width(XLNAD))/1000000
  # nonXLNAD.MB <- sum(width(nonXLNAD))/1000000
  # 
  # NAD.MB <- sum(width(NAD))/1000000
  # nonNAD.MB <- sum(width(nonNAD))/1000000
  
  mbl <- lapply(1:length(XL.nonXL.subset),function(u,XL.nonXL.subset){
    
    mb <- sum(width(XL.nonXL.subset[[u]]))/1000000
    mb
  },XL.nonXL.subset)
  names(mbl) <- names(XL.nonXL.subset)
  #names(mbl) <- unique(as.character(YYYY$SetName))[1:6]
  
  mbl2 <- do.call(rbind,mbl)
  colnames(mbl2) <- "MB"
  
  mbl3 <- as.data.frame(c(XLNAD.MB,XL_nonNAD.MB,nonXLNAD.MB,nonXL_nonNAD.MB,whole.genome.MB))
 
  
  
  row.names(mbl3) <- c("XLNAD.MB","XL_nonNAD.MB","nonXLNAD.MB","nonXL_nonNAD.MB","whole.genome.MB")
  colnames(mbl3) <- "MB"
  mbl4 <- rbind.data.frame(mbl2,mbl3)
  
  YYYYY <- YYYY
    
  YYYY.wholeGenome <- YYYY[which(YYYY$SetName %in% c("wholeGenome")),]
  
  YYYY.XLNAD <- YYYY[which(YYYY$SetName %in% c("XL_Uniq","XL_ciLAD","XL_LAD")),]
  YYYY.XLNAD.3.class <- YYYY.XLNAD
  YYYY.XLNAD$SetName <- "XLNAD" 
    
  YYYY.XL_nonNAD <- YYYY.wholeGenome[-which(YYYY.wholeGenome$GeneName %in% YYYY.XLNAD$GeneName),]
  YYYY.XL_nonNAD$SetName <- "XL_nonNAD" 
  
  YYYY.nonXLNAD <- YYYY[which(YYYY$SetName %in% c("nonXL_Uniq","nonXL_ciLAD","nonXL_LAD")),]
  YYYY.nonXLNAD.3.class <- YYYY.nonXLNAD
  YYYY.nonXLNAD$SetName <- "nonXLNAD" 
  
  YYYY.nonXL_nonNAD <- YYYY.wholeGenome[-which(YYYY.wholeGenome$GeneName %in% YYYY.nonXLNAD$GeneName),]
  YYYY.nonXL_nonNAD$SetName <- "nonXL_nonNAD"
  
  YYYY <- rbind.data.frame(YYYY.XLNAD.3.class,YYYY.XLNAD,YYYY.XL_nonNAD,YYYY.nonXLNAD.3.class,YYYY.nonXLNAD,YYYY.nonXL_nonNAD,YYYY.wholeGenome)
  
  YYYY$SetName <- droplevels(YYYY$SetName)
  
  png(file = file.path(bed.dir,paste0(paste(levels(YYYY$SetName),collapse = "-"),".png")),width = 1500, height = 480)
  boxplot(log10(as.numeric(as.character(YYYY$FPKM))+1)~YYYY$SetName,ylab = "log10(FPKM+1)")
  dev.off()
  
  
  num.gene <- lapply(unique(as.character(YYYY$SetName)), function(u,YYYY){
    num.gene <- length(unique(as.character(YYYY[which(YYYY$SetName==u),]$GeneName)))
    num.gene
  },YYYY)
  names(num.gene) <- unique(as.character(YYYY$SetName))
  num.gene.1 <- do.call(rbind,num.gene)
  colnames(num.gene.1) <- "num"
  
  row.names(mbl4) <- gsub(".MB","",row.names(mbl4))
  row.names(mbl4) <- gsub("whole.genome","wholeGenome",row.names(mbl4))
  
  mbl5 <- merge(mbl4,num.gene.1,by=0)
  
  #XLNAD.num.gene <- length(unique(as.character(YYYY[which(YYYY$SetName %in% c("XL_Uniq","XL_ciLAD","XL_LAD")),]$GeneName)))
  #nonXLNAD.num.gene <- length(unique(as.character(YYYY[which(YYYY$SetName %in% c("nonXL_Uniq","nonXL_ciLAD","nonXL_LAD")),]$GeneName)))
  
  #XLNAD.nonXLNAD.num.gene <- as.data.frame(c(XLNAD.num.gene,nonXLNAD.num.gene))
  #colnames(XLNAD.nonXLNAD.num.gene) <- "num"
  #row.names(XLNAD.nonXLNAD.num.gene) <- c("XLNAD","nonXLNAD")
  
  #num.gene.1 <- rbind(num.gene.1,XLNAD.nonXLNAD.num.gene)
  
  
  #mbl4[rownames(mbl4)=="XLNAD",] <- sum(mbl4[rownames(mbl4) %in% c("XL_Uniq","XL_ciLAD","XL_LAD"),])
  
  #mbl4[rownames(mbl4)=="nonXLNAD",] <- sum(mbl4[rownames(mbl4) %in% c("nonXL_Uniq","nonXL_ciLAD","nonXL_LAD"),])
  
  mb.gene.num <- cbind(mbl5$MB,mbl5$num,mbl5$num/mbl5$MB,mbl5$Row.names)
  
  #mb.gene.num <- data.frame(mb.gene.num,Category=row.names(mb.gene.num))
  
  colnames(mb.gene.num) <- c("MB","num.gene","Density","Category")
  
  DF <- as.data.frame(mb.gene.num)
  DF$MB <- as.numeric(as.character(DF$MB))
  DF$num.gene <- as.numeric(as.character(DF$num.gene))
  DF$Density <- as.numeric(as.character(DF$Density))
  
  #create a colour scheme based on grouping variable 'zone'
  Category<- as.character(DF$Category)
  color.codes <- rainbow(length(Category))
  
  c.1 <- color.codes[which(Category %in% c("nonXL_Uniq","XL_Uniq"))][1]

  color.codes[which(Category %in% c("nonXL_Uniq","XL_Uniq"))] <- c.1
  
  c.2 <- color.codes[which(Category %in% c("nonXL_ciLAD","XL_ciLAD"))][1]
  color.codes[which(Category %in% c("nonXL_ciLAD","XL_ciLAD"))] <- c.2
  
  c.3 <- color.codes[which(Category %in% c("nonXL_LAD" ,"XL_LAD"))][1]
  color.codes[which(Category %in% c("nonXL_LAD" ,"XL_LAD"))] <- c.3
    
  c.4 <- color.codes[which(Category %in% c("nonXLNAD" ,"XLNAD"))][1]
  color.codes[which(Category %in% c("nonXLNAD" ,"XLNAD"))] <- c.4
  
  #color.codes<-as.character(c("#3399FF", "#FF0000","#9633FF","#3399FF", "#FF0000","#9633FF","#0000FF","#C71585","#FFC0CB","#FF0000FF","#00FFFFFF"))
  
  df2=data.frame(Category, color.codes)
  
  # merge color specifications with data
  df <-merge(DF,df2, by=("Category"), all.x=TRUE, all.y=TRUE)
  
  
  df$Category <- factor(df$Category, levels = df$Category[order(df$Density,decreasing = T)])
  
  sp <- ggplot(df, aes(x = Category, y = Density, fill = color.codes)) +
    geom_bar(stat="identity",position = "dodge")  +
    theme(plot.title=element_text(hjust=0.5)) +  
    ylab("Gene density (# of genes/MB)")  +  theme(legend.position="none")
  
  spp <- list(sp=sp)
  multi.page <- ggarrange(plotlist=spp,nrow = 1, ncol = 1)
  
  if(!dir.exists(bed.dir)){dir.create(bed.dir,recursive = TRUE)}
  ggexport(multi.page, width = 1000, height = 480,filename = file.path(bed.dir,"GeneDensity.png"))
  
  write.table(df[,1:4],file = file.path(bed.dir,paste0("Density",".txt")),
              append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T)
  
}

getCTCFDB <- function() {
  LOLA.collection.folder <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/MouseCollection"
  
  db.info <- read.table("/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/MouseCollection/TF_mouse_data_information.txt",header = F,sep ="\t")
  colnames(db.info) <- c("Index","dataSource","Organism","None","cellType","tissue","antibody","filename")
  
  db.info.1 <- db.info[,c(8,5,6,7,2)]
  
  db.info.2 <- db.info[which(db.info.1[,4]=="CTCF"),]
  
  write.table(db.info.1,file = file.path(LOLA.collection.folder,paste0("all",".txt")),
              append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T)
  
  write.table(db.info.2,file = file.path(LOLA.collection.folder,paste0("index",".txt")),
              append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T)
  
  ctcf.files <- file.path(LOLA.collection.folder,"TF_mouse",db.info.2$filename)
  ctcf.Mm.collection <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/MouseCtcf/regions"
  
  if(!dir.exists(ctcf.Mm.collection)){dir.create(ctcf.Mm.collection,recursive = TRUE)}
  
  cmd = paste("cp",ctcf.files,file.path(ctcf.Mm.collection,"."))
  
  null <- lapply(cmd, function(u){system(u)})
  
  
  RegionDB.ctcf <- loadRegionDB("/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin")
  
  output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Results_10_6_with_other_178_CTCF"
  getOverlapWithOther(XL.nonXL.subset,select.query.peak.index= NULL,RegionDB.ctcf,output.file.dir)
  
  output.file.dir <- "/Users/aiminyan/Aimin/ProjectAtCluster/umw_paul_kaufman/Aimin/Results_10_6_2018_ctcf"
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
    
  saveRDS(XL.nonXL.subset,file.path(output.file.dir,"XL_nonXL.rds"))
  saveRDS(RegionDB.ctcf,file.path(output.file.dir,"ctcfRegionDB.rds"))
}

# input.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/BedFiles"
# re.out.rt.mef <- getMoreBedFiles(input.file.dir,re.out)

getMoreBedFiles <- function(input.file.dir,re.out) {
  file.name.bedgraph <- list.files(input.file.dir,pattern="bedgraph",all.files = TRUE,full.names = TRUE,recursive = FALSE,include.dirs = TRUE)
  
  bedgraph.in<-lapply(file.name.bedgraph,function(u){
    
    if(length(grep("Peric",u)) ==1){
      peaks=read.table(u)
    }else if(length(grep("GSM",u)) ==1){
      peaks=read.table(u)
    }else
    {
      peaks=read.table(u,skip=1)
    }
    colnames(peaks)[1:3]= c("chr","start","end")
    peaks=toGRanges(peaks)
    peaks
  })
  
  names(bedgraph.in) <- gsub(" ","_",tools::file_path_sans_ext(basename(file.name.bedgraph)))
  names(bedgraph.in) <- gsub("_mm10_copy","",names(bedgraph.in))
  re.out.rt.mef <- c(re.out,bedgraph.in)
  re.out.rt.mef
}

# Late RT/NAD/LAD
# Early RT/NAD/LAD
# Late RT/NAD/ciLAD
# Early RT/NAD/ciLAD

# output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_10_8_2018_venn_with_RT_MEF"
# Late_RT
# peak.index <- c(10,8,7)
# name <- c("Late_RT","XL","LAD")
# getCount4Venn(re.out.rt.mef,peak.index,name,output.file.dir)

# peak.index <- c(10,6,7)
# name <- c("Late_RT","nonXL","LAD")
# getCount4Venn(re.out.rt.mef,peak.index,name,output.file.dir)

# peak.index <- c(10,8,1)
# name <- c("Late_RT","XL","ciLAD")
# getCount4Venn(re.out.rt.mef,peak.index,name,output.file.dir)

# peak.index <- c(10,6,1)
# name <- c("Late_RT","nonXL","ciLAD")
# getCount4Venn(re.out.rt.mef,peak.index,name,output.file.dir)

# Early_RT
# peak.index <- c(9,8,7)
# name <- c("Early_RT","XL","LAD")
# getCount4Venn(re.out.rt.mef,peak.index,name,output.file.dir)

# peak.index <- c(9,6,7)
# name <- c("Early_RT","nonXL","LAD")
# getCount4Venn(re.out.rt.mef,peak.index,name,output.file.dir)

# peak.index <- c(9,8,1)
# name <- c("Early_RT","XL","ciLAD")
# getCount4Venn(re.out.rt.mef,peak.index,name,output.file.dir)

# peak.index <- c(9,6,1)
# name <- c("Early_RT","nonXL","ciLAD")
# getCount4Venn(re.out.rt.mef,peak.index,name,output.file.dir)


# input.file.dir <- "~/Aimin/DropboxUmass/NADfinder/BedFiles/Aizhan\ F121-9\ comparisons"
# re.out.rt.mef.add.more.bed <-getMoreBedFiles2(input.file.dir,re.out.rt.mef)

getMoreBedFiles2 <- function(input.file.dir,re.out) {
  file.name.bedgraph <- list.files(input.file.dir,pattern="bed",all.files = TRUE,full.names = TRUE,recursive = FALSE,include.dirs = TRUE)
  
  bedgraph.in<-lapply(file.name.bedgraph,function(u){
    
    if(length(grep("Peric",u)) ==1){
      peaks=read.table(u)
    }else if(length(grep("GSM",u)) ==1){
      peaks=read.table(u)
    }else
    {
      peaks=read.table(u,skip=1)
    }
    colnames(peaks)[1:3]= c("chr","start","end")
    peaks=toGRanges(peaks)
    peaks
  })
  
  names(bedgraph.in) <- gsub(" ","_",tools::file_path_sans_ext(basename(file.name.bedgraph)))
  names(bedgraph.in) <- gsub("_mm10_copy","",names(bedgraph.in))
  re.out.rt.mef <- c(re.out,bedgraph.in)
  re.out.rt.mef
}

# ciLAD
# cLAD
# F121-9 ESC NAD
# XL MEF NAD

# output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_10_9_2018_venn_with_F121_9_XL_MEF"

# peak.index <- c(1,7,12,14)
# name <- c("ciLAD","LAD","F121-9","XL_MEF")
# getCount4Venn(re.out.rt.mef.add.more.bed,peak.index,name,output.file.dir)

# peak.index <- c(11,13)
# name <- c("F121-129","F121_cast")
# getCount4Venn(re.out.rt.mef.add.more.bed,peak.index,name,output.file.dir)

# overlap.129.cast <- findOverlapsOfPeaks(re.out.rt.mef.add.more.bed[c(11,13)])
# unique.peaks.in.overlap.129.cast <- getPeaksInUnique(overlap.129.cast)
# output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_10_23_2018_129_cast_unique_peaks_bed"
# outGrl(unique.peaks.in.overlap.129.cast,output.file.dir)



getLADs <- function() {
  GSE17051 <- read.table("/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output/GSE17051_HMM_state_calls_per_probe.txt",header = T)
  
  GSE17051.1 <- GSE17051[which(GSE17051$LaminB1_in_mouse_ESC==1),]
  
  GSE17051.2 <- GSE17051.1[which(GSE17051.1$LaminB1_in_mouse_NPC==0|GSE17051.1$LaminB1_in_mouse_AC==0|GSE17051.1$LaminB1_in_mouse_MEF==0),]
  
  LADs.out.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/Aimin/Output"
  
  write.table(GSE17051.2,file = file.path(LADs.out.dir,paste0("fLADs",".txt")),
              append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T)
  
  GSE17051.3 <- GSE17051[-which(GSE17051$LaminB1_in_mouse_ESC==1),]
  
  GSE17051.4 <- GSE17051.3[which(GSE17051.3$LaminB1_in_mouse_NPC==1|GSE17051.3$LaminB1_in_mouse_AC==1|GSE17051.3$LaminB1_in_mouse_MEF==1),]
  
  write.table(GSE17051.4,file = file.path(LADs.out.dir,paste0("fiLADs",".txt")),
              append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T)
}

# path <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/H3K27me3_bigwig"
# path <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/NR2F2kd_bigwig"

getCvg <- function(path,output.file.dir) {
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  # K27me3.bw.dir <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/H3K27me3_bigwig"
  # NR2F2kd.bw.dir <- "~/Aimin/ProjectAtCluster/umw_nathan_lawson/toDoForAVpaper/1_DensityPlots/NR2F2kd_bigwig"
  
  files <- dir(path, "bw")
  
  if(.Platform$OS.type != "windows"){
    cvglists <- sapply(file.path(path, files), import, 
                       format="BigWig", 
                       as="RleList")
  }else{## rtracklayer can not import bigWig files on Windows
    load(file.path(path, "cvglist.rds"))
  }
  
  saveRDS(cvglists,file = file.path(output.file.dir,"cvglist.rds"))
  
}

# input.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/BedFiles"
# bed.in <- getBedFiles(input.file.dir)

# input.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/BedFiles4GO"
# bed.in <- getBedFiles(input.file.dir)

getBedFiles <- function(input.file.dir) {
  file.name.bedgraph <- list.files(input.file.dir,pattern="*.bed$",all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
  
  bed.in<-lapply(file.name.bedgraph,function(u){
    
    if(length(grep("Peric",u)) ==1){
      peaks=read.table(u)
    }else if(length(grep("GSM",u)) ==1){
      peaks=read.table(u)
    }else
    {
      peaks=read.table(u,skip=1)
    }
    colnames(peaks)[1:3]= c("chr","start","end")
    peaks=toGRanges(peaks)
    peaks
  })
  
  names(bed.in) <- gsub(" ","_",tools::file_path_sans_ext(basename(file.name.bedgraph)))
  names(bed.in) <- gsub("_mm10_copy","",names(bed.in))
  bed.in
  
}

#findOverlapsOfPeaks(re.out.rt.mef.add.more.bed[c(1,2,3)])

#  dir.name="~/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_10_22_2018/nonXL/Bed"
#' input.file.pattern="*.bed$"
#' out.dir.name="~/Aimin/DropboxUmass/NADfinder/Aimin/Output/Results_10_22_2018/nonXL2"
#' txdb="mm10"
#' DD=5000
#'
#' anno.res <- AnntationUsingChipSeeker(dir.name,input.file.pattern,out.dir.name,txdb=txdb,DD,distanceToTSS_cutoff=10000)

AnntationUsingChipSeeker <- function(dir.name,input.file.pattern,out.dir.name,txdb=c("hg19","hg38","mm10"),DD,distanceToTSS_cutoff=5000,assignGenomicAnnotation=TRUE,AP=c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic")) {
  
  file.1 <- list.files(dir.name,pattern=input.file.pattern,all.files = TRUE,full.names = TRUE,recursive = TRUE,include.dirs = TRUE)
  
  re.peaks.only.bed.2 <- file.1
  
  # txdb<-match.arg(txdb)
  
  txdb<-match.arg(txdb,c("hg19","hg38","mm10"),several.ok=TRUE)
  
  
  APpath <- paste(AP,collapse = "_")
  
  temp3=file.path(out.dir.name,"Annotation",APpath)
  
  if(!dir.exists(temp3)){dir.create(temp3,recursive = TRUE)}
  
  d=DD
  res <- lapply(1:length(re.peaks.only.bed.2),function(u,re.peaks.only.bed.2,d){
    u <- 3
    peaks=readPeakFile(re.peaks.only.bed.2[[u]],as="data.frame")
    
    print(head(peaks))
    
    switch (txdb,
            hg38 = {
              cat("Use hg38\n")
              txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
              peakAnno <- annotatePeak(re.peaks.only.bed.2[[u]], tssRegion=c(-d, d),
                                       TxDb=txdb,assignGenomicAnnotation=assignGenomicAnnotation,
                                       genomicAnnotationPriority=AP,annoDb="org.Hs.eg.db")
            },
            hg19 = {
              cat("Use hg19\n")
              txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
              peakAnno <- annotatePeak(re.peaks.only.bed.2[[u]], tssRegion=c(-d, d),
                                       TxDb=txdb,assignGenomicAnnotation=assignGenomicAnnotation,
                                       genomicAnnotationPriority=AP,annoDb="org.Hs.eg.db")
            },
            mm10 = {
              cat("Use mm10\n")
              txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
              peakAnno <- annotatePeak(re.peaks.only.bed.2[[u]], tssRegion=c(-d, d),
                                       TxDb=txdb,assignGenomicAnnotation=assignGenomicAnnotation,
                                       genomicAnnotationPriority=AP,annoDb="org.Mm.eg.db")
            }
            
    )
    
    dropAnnoM <- function (csAnno, distanceToTSS_cutoff)
    {
      idx <- which(abs(mcols(csAnno@anno)[["distanceToTSS"]]) <
                     distanceToTSS_cutoff)
      csAnno@anno <- csAnno@anno[idx]
      csAnno@peakNum <- length(idx)
      if (csAnno@hasGenomicAnnotation) {
        csAnno@annoStat <- ChIPseeker:::getGenomicAnnoStat(csAnno@anno)
        csAnno@detailGenomicAnnotation = csAnno@detailGenomicAnnotation[idx,                                                                        ]
      }
      csAnno
    }
    
    peakAnno <- dropAnnoM(peakAnno,distanceToTSS_cutoff = distanceToTSS_cutoff)
    
    x_name=basename(re.peaks.only.bed.2[[u]])
    cat(x_name)
    
    png(file.path(temp3,paste0(x_name,"_",d,"_around_tss_annotation_pie.png")))
    plotAnnoPie(peakAnno)
    dev.off()
    
    peaks.anno=as.data.frame(peakAnno)
    
    print(head(peaks.anno))
    
    print(colnames(peaks.anno))
    
    peaks.anno.1 <- data.frame(peakName=rownames(peaks.anno),peaks.anno)
    
    write.table(peaks.anno.1,file=file.path(temp3,paste0(x_name,"_",d,"_around_tss_annotation_4_only_mapped_peaks.xls")),
                row.names = FALSE,quote=FALSE,sep="\t")
    
    
    unmapped.peaks<-peaks[-which(paste0(peaks[,1],"_",peaks[,2],"_",peaks[,3]) %in% paste0(peaks.anno[,1],"_",peaks.anno[,2],"_",peaks.anno[,3])),]
    
    cat(dim(peaks)[1]," ",dim(peaks.anno)[1]," ",dim(unmapped.peaks)[1],"\n")
    
    if(dim(unmapped.peaks)[1]!=0){
      
      colnames(unmapped.peaks)=colnames(peaks.anno)[1:3]
      
      unmapped.peaks.3<-smartbind(peaks.anno,unmapped.peaks)
      
      unmapped.peaks.4<-unmapped.peaks.3[order(unmapped.peaks.3[,1],unmapped.peaks.3[,2]),]
      
      unmapped.peaks.4[which(is.na(unmapped.peaks.4$annotation)),6] <- "NotAnno"
      
      row.names(unmapped.peaks.4) <- seq(1,dim(unmapped.peaks.4)[1])
      
      unmapped.peaks.5 <- data.frame(peakName=rownames(unmapped.peaks.4),unmapped.peaks.4)
      
      write.table(unmapped.peaks.5,file=file.path(temp3,paste0(x_name,"_",d,"_around_tss_annotation_4_all_peaks.xls")),row.names = FALSE,quote=FALSE,sep="\t")
      
      f.only <- file.path(temp3,paste0(x_name,"_",d,"_Pie_peaksMappedOnly.pdf"))
      f.all <- file.path(temp3,paste0(x_name,"_",d,"_Pie_peaksAll.pdf"))
      
      makePieChart(peaks.anno,f.only)
      makePieChart(unmapped.peaks.4,f.all)
      
      re <- list(peaksMappedOnly=peaks.anno,peaksAll=unmapped.peaks.4)
      re
    }else
    {
      re <- list(peaksMappedOnly=peaks.anno,peaksAll=peaks.anno)
      re
    }
    
  },re.peaks.only.bed.2,d)
  res
}