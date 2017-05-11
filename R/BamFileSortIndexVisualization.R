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
#' 
#' BamFileSortIndexVisualization(input.file.dir,output.file.dir,genome)
#' 
BamFileSortIndexVisualization <- function(input.file.dir,output.file.dir,genome) {
  
  #library(ChIPpeakAnno)
  
  re<-ParserReadFiles(input.file.dir,"bam",output.file.dir)
  
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

#'R -e 'library(ChipSeq);library(ThreeUTR);ChipSeq:::BamFileSortIndexVisualization2(input.file.dir="/scratch/projects/bbc/Project/Danny_chip2/Alignment/BWA",file.type="marked.bam",output.file.dir="/scratch/projects/bbc/aiminy_project/DannyNewNgsPlot")'

BamFileSortIndexVisualization2 <- function(input.file.dir,file.type,output.file.dir,BigMem = FALSE,cores = 15, Memory = 25000, Wall.time = "72:00", span.ptile = 8) {
  
  #library(ChIPpeakAnno)
  
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
  
  lapply(1:length(re.out),function(u,m.id, Wall.time, cores, Memory,
                                   span.ptile,re.out,temp3){
    
    file_name = file_path_sans_ext(basename(re.out[[u]]))
    
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
      
      job.name = paste0("bamSort.", u)
      cmd1 = paste0("bsub -P bbc -J \"", job.name, paste0("\" -o %J.",
                                                          job.name, ".log "), paste0("-e %J.", job.name, ".err -W"))
      
      
      cmd2=paste("samtools sort",re.out[[u]],file.path(temp3,paste0(file_name, "_sotrted")),sep=" ")
      cmd3 = paste(cmd1, cmd0, cmd2, sep = " ")
    } else
    {
        cmd3 = paste("samtools sort",re.out[[u]],file.path(temp3,paste0(file_name, "_sotrted")),sep=" ")
    }
    
    cmd <- cmd3
    
    system(cmd)
    
    cmd
  
  },m.id,Wall.time,cores,Memory,span.ptile,re.out,temp3)

    
  lapply(1:length(re.out),function(u,m.id, Wall.time, cores, Memory,
                                   span.ptile,re.out,temp3){
    
    file_name = file_path_sans_ext(basename(re.out[[u]]))
    
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
      
      job.name = paste0("bamIndex.", u)
      wait.job.name = paste0("bamSort.", u)
      cmd1 = paste0("bsub -w \"done(\"", wait.job.name, "\")\"", " -P bbc -J \"",
                    job.name, paste0("\" -o %J.", job.name, ".log "), paste0("-e %J.",
                                                                             job.name, ".err -W"))
      cmd2=paste("samtools index",re.out[[u]],file.path(temp3,paste0(file_name, "_sotrted.bam")),sep=" ")
      cmd3 = paste(cmd1, cmd0, cmd2, sep = " ")
    } else
    {
      cmd3 = paste("samtools sort",re.out[[u]],file.path(temp3,paste0(file_name, "_sotrted.bam")),sep=" ")
    }
    
    cmd <- cmd3
    
    system(cmd)
    
    cmd
    
  },m.id,Wall.time,cores,Memory,span.ptile,re.out,temp3)
  
  
  lapply(1:length(re.out),function(u,m.id, Wall.time, cores, Memory,
                                   span.ptile,re.out,temp3){
    
    file_name = file_path_sans_ext(basename(re.out[[u]]))
    
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
      cmd5="ngs.plot.r -G hg19 -R tss -C"
      cmd6="-O"
      cmd7="-L 4000"
      cmd8=paste(cmd5,file.path(temp3,paste0(file_name,"_sorted.bam")),cmd6,file.path(temp3,paste0(file_name,"_sorted")),cmd7,sep=" ")
      
      cmd3 = cmd8
    }
    
    cmd <- cmd3
    
    system(cmd)
    
    cmd
    
  },m.id,Wall.time,cores,Memory,span.ptile,re.out,temp3)
  
}