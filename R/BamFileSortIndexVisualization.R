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