input.rds.dir <- "/Users/aiminyan/Aimin/DropboxUmass/NADfinder/tileCount-NADfinder1.5.2"

file.name <- list.files(input.rds.dir,pattern=".*RDS$",all.files = TRUE,full.names = TRUE,recursive = FALSE,include.dirs = TRUE)

rds.l <- lapply(file.name, function(u){
  se18 <- readRDS(u)
})

names(rds.l) <- basename(file.name)
rds.index <- c(24,26,28)

makeSeList <- function(rds.l,rds.index) {
  
  XL.index <- c(grep(rds.index[1],names(rds.l)),grep(paste0("S_",rds.index[2]),names(rds.l)),grep(rds.index[3],names(rds.l)))
  
  rds.l.1 <- rds.l[XL.index]
  
  #se29 <- readRDS("1.tileCounts.S_29.RDS")
  #se33 <- readRDS("1.tileCounts.S_33.RDS")
  
  se <- do.call(cbind,rds.l[XL.index])
 
  metadata(se)$lib.size.chrom <- cbind(metadata(rds.l.1[[grep(rds.index[1],names(rds.l.1))]])$lib.size.chrom, 
                                       metadata(rds.l.1[[grep(rds.index[2],names(rds.l.1))]])$lib.size.chrom,
                                       metadata(rds.l.1[[grep(rds.index[3],names(rds.l.1))]])$lib.size.chrom)
  
  nuc.cols <- colnames(se)[c(2,4,6)]
  gen.cols <- colnames(se)[c(1,3,5)]
  ## Calculate ratios for nucleoloar vs genomic samples.
  
  se <- log2se(se, 
               nucleolusCols = nuc.cols,
               genomeCols = gen.cols,
               transformation = "log2CPMRatio")
  
  seList <- smoothRatiosByChromosome(se, chrom.level.background = FALSE)
  seList
}


rds.index <- c(24,26,28)
seList.XL <- makeSeList(rds.l,rds.index)

rds.index <- c(18,29,33)
seList <- makeSeList(rds.l,rds.index)

