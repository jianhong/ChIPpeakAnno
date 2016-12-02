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
  
  #library(ChIPpeakAnno)
  
  #dir.name=input.file.dir
  #input.file.pattern=input.file.pattern
  
  #dir.name=reformatPath(dir.name)
  
  #output.dir.name=reformatPath(output.file.dir)
  
  #print(output.dir.name)
  
  #file.name=paste0(dir.name,dir(dir.name,recursive = TRUE,pattern=input.file.pattern))
  
  #file.name.2<-as.list(file.name)
  
  #names(file.name.2)=sapply(strsplit(file.name,split="\\/"),"[[",index.file)
  
  re<-ParserReadFiles(input.file.dir,"bam",output.file.dir)
  
  #print(file.name.2)
  
  file.name.2<-re$input
  output.dir.name=re$output
  
  #temp=Sys.time()
  #temp1=gsub(":","-",temp)
  #temp2=gsub(" ","-",temp1)
  
  temp3=paste0(output.dir.name,"_PeakCall")
  
  dir.create(temp3)
  
  re.out<-file.name.2
  
  # cmd9="macs14 -t "$DIR""$trt_file_name" -f BAM -g hs -n "$results_dir""$trt_file_name"_hs_1.00e-05_macs142 -m 6,18 --bw=200 -B -p 0.00001
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



