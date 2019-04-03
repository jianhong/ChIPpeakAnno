output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/NADSplittingRegion"
if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}

gr.whole.genome <- GRanges(seqinfo(BSgenome.Mmusculus.UCSC.mm10))

other.methods <- names(bed.in)[c(3,6,7,8)]
con.trast.2.out.of.3 <- unique(combn(as.character(other.methods),2,simplify=FALSE))

con.trast.1 <- lapply(con.trast.2.out.of.3, function(u){
  n <- length(unique(u))
  if(n ==1 ){
    x <- NULL
  }
  else{x <- u}
  x
})

print(names(bed.in))

getNADSplittingRegion <- function(gr.whole.genome, bed.in, con.trast.1, which.nad, output.file.dir) {
  
  NAD.peaks.nonXL.rest <- GenomicRanges::setdiff(gr.whole.genome,bed.in[[which.nad]])
  
  null <- lapply(con.trast.1, function(u){
    temp.2.pairs <- bed.in[which(names(bed.in) %in% u)]
    temp.2.pairs.overlapping <- GenomicRanges::intersect(temp.2.pairs[[1]],temp.2.pairs[[2]])
    NADSplittingRegion <- GenomicRanges::intersect(NAD.peaks.nonXL.rest,temp.2.pairs.overlapping)
    output.file.name <- paste(c(names(bed.in[which.nad]),names(temp.2.pairs)),collapse = "_")
    orderPeakAndOutPut(as.data.frame(NADSplittingRegion),output.file.dir,output.file.name,outHeader= FALSE)
  })
  
}

getNADSplittingRegion(gr.whole.genome,bed.in,con.trast.1,4,output.file.dir)
getNADSplittingRegion(gr.whole.genome,bed.in,con.trast.1,5,output.file.dir)

getNADSplittingRegion500 <- function(gr.whole.genome, bed.in, con.trast.1, which.nad, output.file.dir) {
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  nonXL.up.500 <- bed.in[[which.nad]]
  end(nonXL.up.500) <- start(nonXL.up.500) + 500
  
  nonXL.down.500 <- bed.in[[which.nad]]
  start(nonXL.down.500) <- end(nonXL.down.500) - 500
  
  nonXL.NAD.Splitting.500 <- c(nonXL.up.500,nonXL.down.500)
  nonXL.NAD.Splitting.500.sorted <- GenomicRanges::sort(nonXL.NAD.Splitting.500)
  
  NAD.peaks.nonXL.rest <-  nonXL.NAD.Splitting.500.sorted
  
  null <- lapply(con.trast.1, function(u){
    temp.2.pairs <- bed.in[which(names(bed.in) %in% u)]
    temp.2.pairs.overlapping <- GenomicRanges::intersect(temp.2.pairs[[1]],temp.2.pairs[[2]])
    NADSplittingRegion <- GenomicRanges::intersect(NAD.peaks.nonXL.rest,temp.2.pairs.overlapping)
    output.file.name <- paste(c(names(bed.in[which.nad]),names(temp.2.pairs)),collapse = "_")
    orderPeakAndOutPut(as.data.frame(NADSplittingRegion),output.file.dir,output.file.name,outHeader= FALSE)
  })
  
}

output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/NADSplittingRegion/NADSplittingRegion500bp"
getNADSplittingRegion500(gr.whole.genome,bed.in,con.trast.1,4,output.file.dir)
getNADSplittingRegion500(gr.whole.genome,bed.in,con.trast.1,5,output.file.dir)

getNADregions <- function(bed.in,which.nad) {

  chr.label <- unique(seqnames(bed.in[[which.nad]]))

  nad.region.bed <- lapply(1:length(chr.label),function(u,chr.label){
    
    chr.index <- which(seqnames(bed.in[[which.nad]])==chr.label[u])
    chr.bed <- bed.in[[which.nad]][chr.index]
    chr.bed.temp <- chr.bed
    start(chr.bed.temp) <- min(start(chr.bed))
    end(chr.bed.temp) <- max(end(chr.bed))
    chr.bed.temp.1 <- unique(chr.bed.temp)
    chr.bed.temp.1
    
  },chr.label)
  
  nad.region.bed.1 <- do.call(c,nad.region.bed)
  nad.region.bed.1
}

nonXL.nad.region <- getNADregions(bed.in,4)
XL.nad.region <- getNADregions(bed.in,5)

con.trast.2 <- con.trast.1[-grep("MACS",con.trast.1)]

getNADSplittingRegionGreaterThan5k <- function(gr.whole.genome, bed.in, con.trast.1, which.nad, output.file.dir) {
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  print(gr.whole.genome)
  
  NAD.peaks.nonXL.rest <- GenomicRanges::setdiff(gr.whole.genome,bed.in[[which.nad]])
  
  null <- lapply(con.trast.1, function(u){
    temp.2.pairs <- bed.in[which(names(bed.in) %in% u)]
    temp.2.pairs.overlapping <- GenomicRanges::intersect(temp.2.pairs[[1]],temp.2.pairs[[2]])
    NADSplittingRegion <- GenomicRanges::intersect(NAD.peaks.nonXL.rest,temp.2.pairs.overlapping)
    NADSplittingRegionFiltered <- NADSplittingRegion[width(NADSplittingRegion)>50000]
    
    #print(NADSplittingRegionFiltered[width(NADSplittingRegionFiltered)==315001])
    
    output.file.name <- paste(c(names(bed.in[which.nad]),names(temp.2.pairs)),collapse = "_")
    orderPeakAndOutPut(as.data.frame(NADSplittingRegionFiltered),output.file.dir,output.file.name,outHeader= FALSE)
  })
  
}

output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/NADSplittingRegion/NADSplittingRegionFiltered"
getNADSplittingRegionGreaterThan5k(nonXL.nad.region, bed.in, con.trast.2, 4, output.file.dir)
getNADSplittingRegionGreaterThan5k(XL.nad.region, bed.in, con.trast.2, 5, output.file.dir)

# file.name.nad.new <- c("/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/Julie/nonXL_MEF_1.5FC_minp_countF200_qlt0.05 igv.bed","/Users/aiminyan/Aimin/DropboxUmass/pub_old_galaxy/Julie/XL_MEF_1.5FC_minp_countF200_qlt0.05 igv.bed")
# 
# bed.in.new <-lapply(1:length(file.name.nad.new),function(u,file.name.nad.new){
#   if(u %in% c(1,2)){peaks=read.table(file.name.nad.new[u],skip=1)}else{peaks=read.table(file.name.nad.new[u])}
#   
#   colnames(peaks)[1:3]= c("chr","start","end")
#   peaks=toGRanges(peaks)
#   peaks
# },file.name.nad.new)
# 
# names(bed.in.new) <- basename(file.name.nad.new)

getNADSplittingRegionGreaterThan5kChrBased <- function(bed.in, con.trast.1, which.nad,width.filter.value,output.file.dir) {
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
    re <- lapply(con.trast.1, function(u,bed.in,which.nad,width.filter.value){
      
      temp.2.pairs <- bed.in[which(names(bed.in) %in% u)]
      
      chr.label <- unique(as.character(seqnames(bed.in[[which.nad]])))
      
      nad.region.bed <- lapply(1:length(chr.label),function(u,chr.label,bed.in,which.nad,temp.2.pairs,width.filter.value){
        
        chr.index <- which(as.character(seqnames(bed.in[[which.nad]]))==chr.label[u])
        chr.bed <- bed.in[[which.nad]][chr.index]
        chr.bed.temp <- chr.bed
        start(chr.bed.temp) <- min(start(chr.bed))
        end(chr.bed.temp) <- max(end(chr.bed))
        chr.bed.temp.1 <- unique(chr.bed.temp)
      
        NAD.peaks.nonXL.rest <- GenomicRanges::setdiff(chr.bed.temp.1,chr.bed)
      
        chr.index.1 <- which(as.character(seqnames(temp.2.pairs[[1]]))==chr.label[u])
        temp.2.pairs.1.chr <- temp.2.pairs[[1]][chr.index.1]
        
        chr.index.2 <- which(as.character(seqnames(temp.2.pairs[[2]]))==chr.label[u])
        temp.2.pairs.2.chr <- temp.2.pairs[[2]][chr.index.2]
      
        hits1 <- findOverlaps(NAD.peaks.nonXL.rest,temp.2.pairs.1.chr,type=c("within"))
        NAD.peaks.nonXL.rest1 <- NAD.peaks.nonXL.rest[queryHits(hits1)]
        
        hits2 <- findOverlaps(NAD.peaks.nonXL.rest,temp.2.pairs.2.chr,type=c("within"))
        NAD.peaks.nonXL.rest2 <- NAD.peaks.nonXL.rest[queryHits(hits2)]
        
        hits3 <- findOverlaps(NAD.peaks.nonXL.rest,temp.2.pairs.1.chr,type=c("equal"))
        NAD.peaks.nonXL.rest3 <- NAD.peaks.nonXL.rest[queryHits(hits3)]
        
        hits4 <- findOverlaps(NAD.peaks.nonXL.rest,temp.2.pairs.2.chr,type=c("equal"))
        NAD.peaks.nonXL.rest4 <- NAD.peaks.nonXL.rest[queryHits(hits4)]
        
        NAD.peaks.nonXL.rest.sorted <- GenomicRanges::sort(c(NAD.peaks.nonXL.rest1,NAD.peaks.nonXL.rest2,NAD.peaks.nonXL.rest3,NAD.peaks.nonXL.rest4))
        
        temp.2.pairs.overlapping <- GenomicRanges::intersect(temp.2.pairs.1.chr,temp.2.pairs.2.chr)
        NADSplittingRegion <- GenomicRanges::intersect(NAD.peaks.nonXL.rest.sorted,temp.2.pairs.overlapping)
        NADSplittingRegionFiltered <- NADSplittingRegion[width(NADSplittingRegion)>=width.filter.value]
        NADSplittingRegionFiltered
      },chr.label,bed.in,which.nad,temp.2.pairs,width.filter.value)
      
      nad.region.bed.1 <- do.call(c,nad.region.bed)
    
      output.file.name <- paste(c(names(bed.in[which.nad]),names(temp.2.pairs)),collapse = "_")
      orderPeakAndOutPut(as.data.frame(nad.region.bed.1),output.file.dir,output.file.name,outHeader= FALSE)
    
      png(file.path(output.file.dir,paste0(output.file.name,".png")))
      boxplot(width(nad.region.bed.1))
      dev.off()
    
      nad.region.bed.1
      
  },bed.in,which.nad,width.filter.value)
    
    names(re) <- paste0(names(bed.in)[which.nad],"_",unlist(lapply(con.trast.1,function(u){paste(u,collapse = "_")})))
    re
}

output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/NADSplittingRegion/NADSplittingRegionFilteredChrbased2"

getNADSplittingRegionGreaterThan5kChrBased(bed.in,con.trast.2,4,output.file.dir)
getNADSplittingRegionGreaterThan5kChrBased(bed.in,con.trast.2,5,output.file.dir)

output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/NADSplittingRegion/NADSplittingRegionFilteredChrbasedAbove1bp"
re1 <- getNADSplittingRegionGreaterThan5kChrBased(bed.in,con.trast.2,4,1,output.file.dir)
re2 <- getNADSplittingRegionGreaterThan5kChrBased(bed.in,con.trast.2,5,1,output.file.dir)

makeSum <- function(re1,re2,output.file.dir) {
  
  re1.1 <- lapply(1:length(re1), function(u,re1){
    x <- width(re1[[u]])
    y <- data.frame(len=x, name=names(re1)[u])
    y
  },re1)
  
  re1.2 <- do.call(rbind.data.frame,re1.1)
  
  re2.1 <- lapply(1:length(re2), function(u,re2){
    x <- width(re2[[u]])
    y <- data.frame(len=x, name=names(re2)[u])
    y
  },re2)
  
  re2.2 <- do.call(rbind.data.frame,re2.1)
  
  re.11.22 <- rbind(re1.2,re2.2)
  re.11.22$name <- gsub("NADfinder_","",re.11.22$name)
  
  png(file.path(output.file.dir,paste0("NSRs_length_distribution",".png")))
  par(mar=c(14,5,4,2)) 
  boxplot(log10(re.11.22$len)~re.11.22$name,las=2)
  dev.off()
  
  png(file.path(output.file.dir,paste0("NSRs_length_histogram",".png")))
  par(mar=c(14,5,4,2)) 
  hist(log10(re.11.22$len))
  dev.off()
  
  print(summary(re.11.22$len))
  print(summary(log10(re.11.22$len)))
}

con.trast.1.out.of.3 <- unique(combn(as.character(other.methods),1,simplify=FALSE))
con.trast.3 <- con.trast[-grep("MACS",con.trast.1.out.of.3)]

getNADSplittingRegionGreaterThan5kChrBasedOneMethod <- function(bed.in, con.trast.1, which.nad,width.filter.value,output.file.dir) {
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  re <- lapply(con.trast.1, function(u,bed.in,which.nad,width.filter.value){
    
    temp.2.pairs <- bed.in[which(names(bed.in) %in% u)]
    
    chr.label <- unique(as.character(seqnames(bed.in[[which.nad]])))
    
    nad.region.bed <- lapply(1:length(chr.label),function(u,chr.label,bed.in,which.nad,temp.2.pairs,width.filter.value){
      
      chr.index <- which(as.character(seqnames(bed.in[[which.nad]]))==chr.label[u])
      chr.bed <- bed.in[[which.nad]][chr.index]
      chr.bed.temp <- chr.bed
      start(chr.bed.temp) <- min(start(chr.bed))
      end(chr.bed.temp) <- max(end(chr.bed))
      chr.bed.temp.1 <- unique(chr.bed.temp)
      
      NAD.peaks.nonXL.rest <- GenomicRanges::setdiff(chr.bed.temp.1,chr.bed)
      
      chr.index.1 <- which(as.character(seqnames(temp.2.pairs[[1]]))==chr.label[u])
      temp.2.pairs.1.chr <- temp.2.pairs[[1]][chr.index.1]
      
    #  chr.index.2 <- which(as.character(seqnames(temp.2.pairs[[2]]))==chr.label[u])
    #  temp.2.pairs.2.chr <- temp.2.pairs[[2]][chr.index.2]
      
      hits1 <- findOverlaps(NAD.peaks.nonXL.rest,temp.2.pairs.1.chr,type=c("within"))
      NAD.peaks.nonXL.rest1 <- NAD.peaks.nonXL.rest[queryHits(hits1)]
      
     # hits2 <- findOverlaps(NAD.peaks.nonXL.rest,temp.2.pairs.2.chr,type=c("within"))
      #NAD.peaks.nonXL.rest2 <- NAD.peaks.nonXL.rest[queryHits(hits2)]
      
      hits3 <- findOverlaps(NAD.peaks.nonXL.rest,temp.2.pairs.1.chr,type=c("equal"))
      NAD.peaks.nonXL.rest3 <- NAD.peaks.nonXL.rest[queryHits(hits3)]
      
      #hits4 <- findOverlaps(NAD.peaks.nonXL.rest,temp.2.pairs.2.chr,type=c("equal"))
      #NAD.peaks.nonXL.rest4 <- NAD.peaks.nonXL.rest[queryHits(hits4)]
      
      #NAD.peaks.nonXL.rest.sorted <- GenomicRanges::sort(c(NAD.peaks.nonXL.rest1,NAD.peaks.nonXL.rest2,NAD.peaks.nonXL.rest3,NAD.peaks.nonXL.rest4))
      
      NAD.peaks.nonXL.rest.sorted <- GenomicRanges::sort(c(NAD.peaks.nonXL.rest1,NAD.peaks.nonXL.rest3))
      
    #  temp.2.pairs.overlapping <- GenomicRanges::intersect(temp.2.pairs.1.chr,temp.2.pairs.2.chr)
      
      temp.2.pairs.overlapping <- temp.2.pairs.1.chr
      
      NADSplittingRegion <- GenomicRanges::intersect(NAD.peaks.nonXL.rest.sorted,temp.2.pairs.overlapping)
      NADSplittingRegionFiltered <- NADSplittingRegion[width(NADSplittingRegion)>=width.filter.value]
      NADSplittingRegionFiltered
    },chr.label,bed.in,which.nad,temp.2.pairs,width.filter.value)
    
    nad.region.bed.1 <- do.call(c,nad.region.bed)
    
    output.file.name <- paste(c(names(bed.in[which.nad]),names(temp.2.pairs)),collapse = "_")
    orderPeakAndOutPut(as.data.frame(nad.region.bed.1),output.file.dir,output.file.name,outHeader= FALSE)
    
    png(file.path(output.file.dir,paste0(output.file.name,".png")))
    boxplot(width(nad.region.bed.1))
    dev.off()
    
    nad.region.bed.1
    
  },bed.in,which.nad,width.filter.value)
  
  names(re) <- paste0(names(bed.in)[which.nad],"_",unlist(lapply(con.trast.1,function(u){paste(u,collapse = "_")})))
  re
}

output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/NADSplittingRegion/NADSplittingRegionFilteredChrbasedAbove1bpOneMethod"
re1.one <- getNADSplittingRegionGreaterThan5kChrBasedOneMethod(bed.in,con.trast.3,4,1,output.file.dir)
re2.one <- getNADSplittingRegionGreaterThan5kChrBasedOneMethod(bed.in,con.trast.3,5,1,output.file.dir)

makeSum(re1.one,re2.one,output.file.dir)

NAD.peaks.nonXL.in.L <- lapply(bed.in[c(3,6,8)],function(u,bed.in){
  
  hits1 <- findOverlaps(bed.in[[4]],u,type=c("within"))
  NAD.peaks.nonXL.within <- bed.in[[4]][queryHits(hits1)]
  
  hits3 <- findOverlaps(bed.in[[4]],u,type=c("equal"))
  NAD.peaks.nonXL.equal <- bed.in[[4]][queryHits(hits3)]
  
  NAD.peaks.nonXL.in <- c(NAD.peaks.nonXL.within,NAD.peaks.nonXL.equal)
  
  NAD.peaks.nonXL.in
  
},bed.in)

do.call(c,NAD.peaks.nonXL.in.L)

c(NAD.peaks.nonXL.in.L[[1]],NAD.peaks.nonXL.in.L[[2]],NAD.peaks.nonXL.in.L[[3]])

getNADSplittingRegionGreaterThan5kChrBased3 <- function(bed.in, con.trast.1, which.nad,width.filter.value,output.file.dir) {
  
  if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}
  
  re <- lapply(con.trast.1, function(u,bed.in,which.nad,width.filter.value){
    
    temp.2.pairs <- bed.in[which(names(bed.in) %in% u)]
    
    #temp.2.pairs[[1]] <- GenomicRanges::reduce(temp.2.pairs[[1]])
    #temp.2.pairs[[2]] <- GenomicRanges::reduce(temp.2.pairs[[2]])
    
    chr.label <- unique(as.character(seqnames(bed.in[[which.nad]])))
    
    nad.region.bed <- lapply(1:length(chr.label),function(u,chr.label,bed.in,which.nad,temp.2.pairs,width.filter.value){

      hits1 <- findOverlaps(bed.in[[which.nad]],temp.2.pairs[[1]],type=c("within"))
      NAD.peaks.nonXL.1 <- bed.in[[which.nad]][queryHits(hits1)]
      
      hits2 <- findOverlaps(bed.in[[which.nad]],temp.2.pairs[[2]],type=c("within"))
      NAD.peaks.nonXL.2 <- bed.in[[which.nad]][queryHits(hits2)]
      
      hits3 <- findOverlaps(bed.in[[which.nad]],temp.2.pairs[[1]],type=c("equal"))
      NAD.peaks.nonXL.3 <- bed.in[[which.nad]][queryHits(hits3)]
      
      hits4 <- findOverlaps(bed.in[[which.nad]],temp.2.pairs[[2]],type=c("equal"))
      NAD.peaks.nonXL.4 <- bed.in[[which.nad]][queryHits(hits4)]
      
      NAD.peaks.nonXL.in <- unique(GenomicRanges::sort(c(NAD.peaks.nonXL.1,NAD.peaks.nonXL.2,NAD.peaks.nonXL.3,NAD.peaks.nonXL.4)))
      
      #print(NAD.peaks.nonXL.in[as.character(seqnames(NAD.peaks.nonXL.in))=="chr19"])
      
      chr.index <- which(as.character(seqnames(NAD.peaks.nonXL.in))==chr.label[u])
      
      if(length(chr.index)!=0){
      chr.bed <- NAD.peaks.nonXL.in[chr.index]
       
      #print(chr.bed)
      chr.bed.temp <- chr.bed
      start(chr.bed.temp) <- min(start(chr.bed))
      end(chr.bed.temp) <- max(end(chr.bed))
      chr.bed.temp.1 <- unique(chr.bed.temp)
    
      NAD.peaks.nonXL.rest <- GenomicRanges::setdiff(chr.bed.temp.1,chr.bed)
      
      chr.index.1 <- which(as.character(seqnames(temp.2.pairs[[1]]))==chr.label[u])
      temp.2.pairs.1.chr <- temp.2.pairs[[1]][chr.index.1]
      
      chr.index.2 <- which(as.character(seqnames(temp.2.pairs[[2]]))==chr.label[u])
      temp.2.pairs.2.chr <- temp.2.pairs[[2]][chr.index.2]
      
      hits1 <- findOverlaps(NAD.peaks.nonXL.rest,temp.2.pairs.1.chr,type=c("within"))
      NAD.peaks.nonXL.rest1 <- NAD.peaks.nonXL.rest[queryHits(hits1)]
      
      hits2 <- findOverlaps(NAD.peaks.nonXL.rest,temp.2.pairs.2.chr,type=c("within"))
      NAD.peaks.nonXL.rest2 <- NAD.peaks.nonXL.rest[queryHits(hits2)]
      
      hits3 <- findOverlaps(NAD.peaks.nonXL.rest,temp.2.pairs.1.chr,type=c("equal"))
      NAD.peaks.nonXL.rest3 <- NAD.peaks.nonXL.rest[queryHits(hits3)]
      
      hits4 <- findOverlaps(NAD.peaks.nonXL.rest,temp.2.pairs.2.chr,type=c("equal"))
      NAD.peaks.nonXL.rest4 <- NAD.peaks.nonXL.rest[queryHits(hits4)]
      
      NAD.peaks.nonXL.rest.sorted <- GenomicRanges::sort(c(NAD.peaks.nonXL.rest1,NAD.peaks.nonXL.rest2,NAD.peaks.nonXL.rest3,NAD.peaks.nonXL.rest4))
      
      temp.2.pairs.overlapping <- GenomicRanges::intersect(temp.2.pairs.1.chr,temp.2.pairs.2.chr)
      NADSplittingRegion <- GenomicRanges::intersect(NAD.peaks.nonXL.rest.sorted,temp.2.pairs.overlapping)
      NADSplittingRegionFiltered <- NADSplittingRegion[width(NADSplittingRegion)>=width.filter.value]
      NADSplittingRegionFiltered
      
     # print(NADSplittingRegionFiltered)
      
      }
    },chr.label,bed.in,which.nad,temp.2.pairs,width.filter.value)
    
    cat(names(temp.2.pairs),"\n")
    # print(nad.region.bed)
    
    cat(length(nad.region.bed),"\n")
    print(nad.region.bed)
    
    # print(which(sapply(nad.region.bed, is.null)))
    # 
    # nad.region.bed.net <- nad.region.bed[-which(sapply(nad.region.bed, is.null))]
    # 
    # print(nad.region.bed.net)
    
    
#    nad.region.bed.net <- lapply(nad.region.bed,function(u)
#      {
#       if(length(u)!=0)
#       {x <- u}
#      x})
     
    #if(length(nad.region.bed)!=0){
    nad.region.bed.1 <- do.call(c,nad.region.bed)
    
    #cat(names(temp.2.pairs),"\n")
    cat(length(nad.region.bed.1),"\n")
    print(nad.region.bed.1)
    
    
    
    #}else{
    #  print(names(temp.2.pairs))
    #}
    output.file.name <- paste(c(names(bed.in[which.nad]),names(temp.2.pairs)),collapse = "_")
    
    if(length(nad.region.bed.1) > 0){
      cat(output.file.name," have ",length(nad.region.bed.1),"NSRs \n\n")
    #  print(length(nad.region.bed.1))
    orderPeakAndOutPut(as.data.frame(nad.region.bed.1),output.file.dir,output.file.name,outHeader= FALSE)
    
    #cat(output.file.name,"\n")
    #print(nad.region.bed.1)
    
    #if(length(nad.region.bed.1) > 0){
    
    png(file.path(output.file.dir,paste0(output.file.name,".png")))
    boxplot(width(nad.region.bed.1))
    dev.off()
    re <- nad.region.bed.1
    }else
    {
      cat(output.file.name," does not have NSRs \n")
    #  print(nad.region.bed.1)
      re <- NULL
    }
  #  cat(output.file.name,"\n")
  #  print(nad.region.bed.1)
    
    #}
    re
    
  },bed.in,which.nad,width.filter.value)
  
  names(re) <- paste0(names(bed.in)[which.nad],"_",unlist(lapply(con.trast.1,function(u){paste(u,collapse = "_")})))
  re
}

output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/NADSplittingRegion/NADSplittingRegionFilteredChrbasedAbove1bpUseNestedNAD"

sink("~/Aimin/DropboxUmass/NADfinder/Aimin/Output/NADSplittingRegion/NADSplittingRegionFilteredChrbasedAbove1bpUseNestedNAD/test.txt")
re1 <- getNADSplittingRegionGreaterThan5kChrBased3(bed.in,con.trast.2[1],4,1,output.file.dir)
sink()
re2 <- getNADSplittingRegionGreaterThan5kChrBased3(bed.in,con.trast.2,5,1,output.file.dir)

re1[[1]] <- NULL
re2[[1]] <- NULL
makeSum(re1,re2,output.file.dir)

nonNAD.nonXL <- GenomicRanges::setdiff(gr.whole.genome,bed.in[[4]])
nonNAD.XL <- GenomicRanges::setdiff(gr.whole.genome,bed.in[[5]])

NSRs <- c(re1,re2)
NSRs.anno <- getAnnotatedGene(NSRs,"Mm")

fpkm.value.ref <- data.frame(SetName=rep("wholeGenome",dim(fpkm.value)[1]),GeneName=fpkm.value$gene.name,FPKM=fpkm.value$fpkm)

getFPKM4NAD <- function(bed.in,which.nad,nad.name) {
  YYY <- getAnnotatedGene(bed.in[which.nad],"Mm")
  NAD.fpkm.value <- fpkm.value.ref[which(fpkm.value.ref$GeneName %in% YYY[[1]]$gene_name),]
  NAD.fpkm.value$SetName <- nad.name
  nonNAD.fpkm.value <- fpkm.value.ref[-which(fpkm.value.ref$GeneName %in% YYY[[1]]$gene_name),]
  nonNAD.fpkm.value$SetName <- paste0("nonNAD.",nad.name)
  YYYY <- rbind(NAD.fpkm.value,nonNAD.fpkm.value)
  YYYY
}

nonXL.fpkm <- getFPKM4NAD(bed.in,4,"NAD-nonXL")
XL.fpkm <- getFPKM4NAD(bed.in,5,"NAD-XL")

NSRs.anno <- getAnnotatedGene(NSRs,"Mm")

NSRs.fpkm <- getFPKM4DiffSet(NSRs.anno,fpkm.value)

NSRs.fpkm <- rbind(NSRs.fpkm,nonXL.fpkm,XL.fpkm,fpkm.value.ref)

type.1.2.nad.anno <- getAnnotatedGene(bed.in.new.4.JC,"Mm")

type.1.2.nad.anno.fpkm <- getFPKM4DiffSet(type.1.2.nad.anno,fpkm.value)

NSRs.fpkm <- rbind(NSRs.fpkm,type.1.2.nad.anno.fpkm)

NSRs.fpkm$FPKM <- log10(as.numeric(as.character(NSRs.fpkm$FPKM))+1)

saveRDS(NSRs.fpkm,file=file.path(output.file.dir,"NSRs_FPKM.rds"))

png(file.path(output.file.dir,paste0("NSRs_FPKM_distribution",".png")))
par(mar=c(18,5,4,2)) 
boxplot(NSRs.fpkm$FPKM~NSRs.fpkm$SetName,las=2)
dev.off()

output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/NADSplittingRegion/NADSplittingRegionFilteredChrbasedAbove1bpUseNestedNAD2"

re1 <- getNADSplittingRegionGreaterThan5kChrBased3(bed.in,con.trast.2,4,1,output.file.dir)
re2 <- getNADSplittingRegionGreaterThan5kChrBased3(bed.in,con.trast.2,5,1,output.file.dir)

output.file.dir <- "~/Aimin/DropboxUmass/NADfinder/Aimin/Output/NADSplittingRegion/NADSplittingRegionFilteredChrbasedAbove1bpUseNestedNADRefinedTest"

sink("~/Aimin/DropboxUmass/NADfinder/Aimin/Output/NADSplittingRegion/NADSplittingRegionFilteredChrbasedAbove1bpUseNestedNADRefinedTest/test_output.txt")
re1 <- getNADSplittingRegionGreaterThan5kChrBased3(bed.in.refined,con.trast.2,4,1,output.file.dir)
re2 <- getNADSplittingRegionGreaterThan5kChrBased3(bed.in.refined,con.trast.2,5,1,output.file.dir)
sink()

makeSum(re1,re2,output.file.dir)

NSRs <- c(re1,re2)
NSRs.anno <- getAnnotatedGene(NSRs,"Mm")

fpkm.value.ref <- data.frame(SetName=rep("wholeGenome",dim(fpkm.value)[1]),GeneName=fpkm.value$gene.name,FPKM=fpkm.value$fpkm)

getFPKM4NAD <- function(bed.in,which.nad,nad.name) {
  YYY <- getAnnotatedGene(bed.in[which.nad],"Mm")
  NAD.fpkm.value <- fpkm.value.ref[which(fpkm.value.ref$GeneName %in% YYY[[1]]$gene_name),]
  NAD.fpkm.value$SetName <- nad.name
  nonNAD.fpkm.value <- fpkm.value.ref[-which(fpkm.value.ref$GeneName %in% YYY[[1]]$gene_name),]
  nonNAD.fpkm.value$SetName <- paste0("nonNAD.",nad.name)
  YYYY <- rbind(NAD.fpkm.value,nonNAD.fpkm.value)
  YYYY
}

nonXL.fpkm <- getFPKM4NAD(bed.in.refined,4,"NAD-nonXL")
XL.fpkm <- getFPKM4NAD(bed.in.refined,5,"NAD-XL")

#NSRs.anno <- getAnnotatedGene(NSRs,"Mm")

NSRs.fpkm <- getFPKM4DiffSet(NSRs.anno,fpkm.value)

NSRs.fpkm <- rbind(NSRs.fpkm,nonXL.fpkm,XL.fpkm,fpkm.value.ref)

type.1.2.nad.anno <- getAnnotatedGene(bed.in.new.4.JC,"Mm")

type.1.2.nad.anno.fpkm <- getFPKM4DiffSet(type.1.2.nad.anno,fpkm.value)

NSRs.fpkm <- rbind(NSRs.fpkm,type.1.2.nad.anno.fpkm)

NSRs.fpkm$FPKM <- log10(as.numeric(as.character(NSRs.fpkm$FPKM))+1)

saveRDS(NSRs.fpkm,file=file.path(output.file.dir,"NSRs_FPKM.rds"))

png(file.path(output.file.dir,paste0("NSRs_FPKM_distribution",".png")))
par(mar=c(18,5,4,2)) 
boxplot(NSRs.fpkm$FPKM~NSRs.fpkm$SetName,las=2)
dev.off()

getTtestP4FPKM <- function(fpkm1,fpkm2){
 x <-  t.test(fpkm1,fpkm2)$p.value
 x  
}

t.p <- sapply(unique(as.character(NSRs.fpkm$SetName)), function(x) sapply(unique(as.character(NSRs.fpkm$SetName)), function(y) getTtestP4FPKM(fpkm1=NSRs.fpkm[as.character(NSRs.fpkm$SetName)==y,]$FPKM, fpkm2=NSRs.fpkm[as.character(NSRs.fpkm$SetName)==x,]$FPKM)))

write.table(t.p,file = file.path(output.file.dir,paste0("FPKM_comparision","_p_value_t_test",".txt")),append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names=NA)

save.image(file = file.path(output.file.dir,"NADSplittingRegion.RData"))
savehistory(file = file.path(output.file.dir,"NADSplittingRegion.Rhistory"))

save.image(file = file.path(output.file.dir,"NADSplittingRegionFPKM.RData"))
savehistory(file = file.path(output.file.dir,"NADSplittingRegionFPKM.Rhistory"))

save.image(file = file.path(output.file.dir,"NADSplittingRegionFPKMWithTypeIAndII.RData"))
savehistory(file = file.path(output.file.dir,"NADSplittingRegionFPKMWithTypeIAndII.Rhistory"))

save.image(file = file.path(output.file.dir,"NADSplittingRegionFPKMWithTypeIAndIIRefined.RData"))
savehistory(file = file.path(output.file.dir,"NADSplittingRegionFPKMWithTypeIAndIIRefined.Rhistory"))

save.image(file = file.path(output.file.dir,"NADSplittingRegionFPKMWithTypeIAndIIRefinedTest.RData"))
savehistory(file = file.path(output.file.dir,"NADSplittingRegionFPKMWithTypeIAndIIRefinedTest.Rhistory"))

save.image(file = file.path(output.file.dir,"NADSplittingRegionFPKMWithTypeIAndIIRefinedTest-4-3-2019.RData"))
savehistory(file = file.path(output.file.dir,"NADSplittingRegionFPKMWithTypeIAndIIRefinedTest-4-3-2019.Rhistory"))
