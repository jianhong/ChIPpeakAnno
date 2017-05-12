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
    
    cellType <- unique(as.character(re11[re11$Type_TF == "Input", 
        ]$Type_Cell))
    
    
    yy <- lapply(cellType, function(u, re11)
    {
        n <- nchar(u)
        index1 <- which(substr(re11$Type_Cell, 1, n) == u)
        z <- re11[index1, ]
        index2 <- which(nchar(as.character(z$Type_Cell)) <= n + 
            2)
        
        z2 <- z[index2, ]
        
        cc <- chipseqCell(u, z2)
        
        cc
    }, re11)
    
    colnames(re1)[4] = "Cell_TF"
    
    Cell.Ab.type <- unique(as.character(re1[, 4]))
    
    re2 <- lapply(1:length(Cell.Ab.type), function(u, Cell.Ab.type, 
        re1)
        {
        
        x = Cell.Ab.type[u]
        z = re1
        
        ZZ <- as.character(z[which(z[, 4] == x), 1])
        ZZ
        
    }, Cell.Ab.type, re1)
    
    names(re2) <- Cell.Ab.type
    
    re21 <- re2[lapply(re2, length) == 1]
    
    re22 <- re2[lapply(re2, length) == 2]
    
    re3 <- list(re1 = re1, re2 = re2, re21 = re21, re22 = re22, 
        re11 = re11, y = yy)
    
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
peakcallwithinput <- function(input.sample.file, input.bam.file, 
    genome = c("Hs", "hs", "HS", "hS"), output.dir, peakcaller = c("macs14", 
        "macs2"), peakPvalue)
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
    
    cellInfo.run <- lapply(1:length(cellInfo), function(u, cellInfo, 
        temp3)
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
            xx.name = paste(ID, gsub(" ", "-", Type_Cell), Type_TF, 
                sep = "-")
            
            cmd12 = paste(cmd9, xx, "-c", x.input, cmd10, file.path(temp3, 
                paste0(xx.name, "_hs_1.00e-05_", peakcaller)), 
                cmd11, sep = " ")
            
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

select.sample <- c("MDA MB 231-DD-1_cJun", "MDA MB 231-1_cJun", 
    "1833-1_cJun")
output.config.dir <- "~/"

configAndMultiplot <- function(res, select.sample, output.config.dir)
{
    
    if (!dir.exists(output.config.dir))
    {
        dir.create(output.config.dir)
    }
    
    x <- res$re11
    
    xx <- x[which(x$Cell_TF %in% select.sample), ]
    
    xxx <- cbind.data.frame(xx$file.name, rep("-1", dim(xx)[1]), 
        gsub(" ", "-", xx$Cell_TF))
    
    config.sample.name <- paste(gsub(" ", "-", select.sample), 
        collapse = "-and-")
    
    config.file.name <- file.path(output.config.dir, paste0(config.sample.name, 
        "-config.txt"))
    
    write.table(xxx, file = config.file.name, col.names = FALSE, 
        row.names = FALSE, quote = FALSE, sep = "\t")
    
    cmd0 <- "ngs.plot.r -G hg19 -R tss -C"
    
    cmd1 <- "-O"
    
    cmd2 <- "-L 4000 -RR 1 -CD 1 -CO \"blue\""
    
    output.results <- file.path(output.config.dir, paste0(config.sample.name, 
        "results"))
    
    cmd <- paste(cmd0, config.file.name, cmd1, output.results, 
        cmd2, sep = " ")
    
    system(cmd)
    
    return(config.file.name)
    
}

#' 
#' input.sample.file <- '/Volumes/Bioinformatics$/2017/DannyNewData/SampleID_INFO_ChIP_new_Danny.csv'
#' 
#' input.bam.file <- '/Volumes/Bioinformatics$/2017/DannyNewData/NewRe2Danny/sorted_bam_files_2.txt'
#' 
#' re <- ChipSeq:::matchBamInputGene(input.sample.file,input.bam.file,"$HOME/all_common_gene_unique.txt","$HOME/NgsConfigFile")
#' 
#' re <- ChipSeq:::matchBamInputGene("/scratch/projects/bbc/aiminy_project/DannyNewData2/SampleID_INFO_ChIP_new_Danny.csv","/scratch/projects/bbc/aiminy_project/DannyNewData2/sorted_bam_files_2.txt","$HOME/all_common_gene_unique.txt","$HOME/NgsConfigFile")
#' 
#' 
matchBamInputGene <- function(input.sample.file, input.bam.file,input.gene.list,output.dir,ngs.para=c("hg19",4000,1,1,"total"))
{
  
  re <- GetSampleInfo(input.sample.file, input.bam.file)
  
  cellInfo <- re$y
  
  output.dir.name = dirname(input.sample.file)
  
  if (!dir.exists(output.dir))
  {
    dir.create(output.dir,recursive = TRUE)
  }
  
  temp3 = output.dir
  
    cmd9 = "ngs.plot.r -G"
    cmd10 = "-R"
    cmd11 = "-C"
    cmd12 = "-O" 
    cmd13 = "-T"
    cmd14 = "-L" 
    cmd15 = "-RR"
    cmd16 = "-CD"
    cmd17= "-GO"
    

  
  cellInfo.run <- lapply(1:length(cellInfo), function(u, cellInfo, 
                                                      temp3)
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
      xx.name = paste(ID, gsub(" ", "-", Type_Cell), Type_TF, 
                      sep = "-")
      
      cmd12 = paste(paste0(xx, ":", x.input),input.gene.list,paste0(gsub(" ","-",Cell_TF), "_cJunAndp27"), sep = "\t")
      cmd12
      cat(cmd12, file = file.path(temp3,paste0(xx.name,"_config_cJunAndp27.txt")),sep = "\t")
    })
    
    x.run
    
  }, cellInfo, temp3)
  
  
  dir.name=temp3
  dir.name=reformatPath(dir.name)
  
  file.name=file.path(dir.name,dir(dir.name,recursive = TRUE))
  
  file.name.2<-as.list(file.name)
  

  # names(file.name.2) = unlist(lapply(file.name.2, function(u)
  # {
  #   u$name
  # }))
  
  zzz <- unlist(file.name.2)
  
  lapply(1:length(zzz), function(u, zzz)
  {
    
    dir.name = dirname(zzz[u][[1]])
    file_name = file_path_sans_ext(basename(zzz[u][[1]]))
    
    cmd = paste("ngs.plot.r -G hg19 -R tss -C",zzz[u][[1]],"-O",file.path(dir.name,paste0(file_name,".tss")),"-T",file_name,"-L 4000 -RR 1 -CD 1 -GO total",sep=" ")
  
    
    
    #system(as.character(zzz[u][[1]]))
    job.name = paste0("bamPlot.", u)
    cmd.pegasus = usePegasus("general", Wall.time = "72:00",cores = 32,Memory = 16000,span.ptile = 16,job.name)
    cmd2 = paste(cmd.pegasus,cmd,sep = " ")
    
    cat(cmd2, "\n")
    cat("\n")
    system(cmd2)
  }, zzz)
  
  
  re <- list(cellInforun = cellInfo.run, zzz = zzz)
  
  #AnntationUsingChipSeeker(temp3, "peaks.bed", temp3, DD = 5000)
  
  return(re)
  
}

#' x <- ChipSeq:::usePegasus("general", Wall.time = "72:00",cores = 32,Memory = 16000,span.ptile = 16,"bamPlot")

usePegasus <- function(job.option=c("general","parallel","bigmem")
                       , Wall.time, cores, Memory, span.ptile,job.name,wait.job.name=NULL) {
  
  job.option <- match.arg(job.option)
  
  switch (job.option,
          parallel = {
            cmd0 = paste(Wall.time,"-n",cores,"-q parallel -R 'rusage[mem=",Memory,"] span[ptile=",span.ptile,"]' -u aimin.yan@med.miami.edu",sep = " ")        
          },
          bigmem = {
            cmd0 = paste(Wall.time,"-n",cores,"-q bigmem -R 'rusage[mem=",Memory,"] span[ptile=", span.ptile, "]' -u aimin.yan@med.miami.edu",sep = " ")
          },
          general = {
            cmd0 = paste(Wall.time,"-n",cores,"-q general -R 'rusage[mem=",Memory,"] span[ptile=",span.ptile, "]' -u aimin.yan@med.miami.edu",sep = " ")
          }
  )
  
  if(!is.null(wait.job.name)){
    cmd1 = paste0("bsub -w \"done(\"", wait.job.name, "\")\"", " -P bbc -J \"",
                  job.name, paste0("\" -o %J.", job.name, ".log "), paste0("-e %J.",
                                                                           job.name, ".err -W"))
  }else{
    cmd1 = paste0("bsub -P bbc -J \"",job.name, paste0("\" -o %J.", job.name, ".log "), paste0("-e %J.",
                                                                                               job.name, ".err -W"))
  }
  
    cmd = paste(cmd1,cmd0,sep=" ")
  
  return(cmd)
}