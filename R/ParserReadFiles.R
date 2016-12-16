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

  file.name=paste0(dir.name,dir(dir.name,recursive = TRUE))
  
  file.name.2<-as.list(file.name)

  file.name.3<-lapply(1:length(file.name.2),function(u,input.file.type,file.name.2){
    
    tmp=file.name.2
    
    x=tmp[[u]]
    path_name=dirname(x)
    file_name=basename(x)
    
    if(file_ext(file_name)==input.file.type){
    re<-paste0(path_name,"/",file_name)
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
  output.dir.name=reformatPath(output.file.dir)
  #print(output.dir.name)
  
  temp=Sys.time()
  temp1=gsub(":","-",temp)
  temp2=gsub(" ","-",temp1)
  temp3=paste0(output.dir.name,"ReadBam_at_",temp2)
  #dir.create(temp3)
  
  re2<-list(input=file.name.4,output=temp3)
  
  return(re2)
}
