#!/usr/bin/env Rscript

library(stringr)

input.csv.file <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/Feston/chromium-dna-sample-indexes-plate.csv"
output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/nathan_lawson/Feston/Output"

if(!dir.exists(output.file.dir)){dir.create(output.file.dir,recursive = TRUE)}

input.csv.data <- read.csv(input.csv.file,header = F)

Y <- lapply(1:dim(input.csv.data)[1],function(u,input.csv.data){
  
  x <- input.csv.data[u,]
  index <- t(x[1])
  index1 <- str_sub(index,7,7)
  index2 <- formatC(as.integer(str_sub(index,8,str_length(index))),width=2, flag="0")
  index3 <- paste0(index1,index2)
  sequence=t(x[-1])
  index4=paste0(index3,"-",formatC(seq(1,dim(sequence)[1]),width=2,flag="0"))
  y <- data.frame(index=index4,sequence=sequence)
  colnames(y) <- c("index","sequence")
  y
  
},input.csv.data)

YY <- do.call(rbind,Y)

write.table(YY,file = file.path(output.file.dir,paste0("reformated",".csv")),append = FALSE, quote = F, sep = ",",eol = "\n", na = "NA", dec = ".", row.names = F)

save.image(file = file.path(output.file.dir,"Feston_reforamt.RData"))
savehistory(file = file.path(output.file.dir,"Feston_reformat.Rhistory"))
