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
#' testR<-ParserReadFiles(input.file.dir,out.file.dir) 
#' 
#' 
ParserReadFiles <- function(input.file.dir,input.file.type,output.file.dir) {
  
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