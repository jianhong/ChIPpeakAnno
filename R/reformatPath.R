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
    
    temp=unlist(strsplit(dir.name,split="\\/"))
    temp[2]="Volumes"
    temp[3]="Bioinformatics$"
    dir.name=paste0(paste0(temp,collapse = "/"),"/")
  }
  
  return(dir.name)
}