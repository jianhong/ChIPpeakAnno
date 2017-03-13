#' useHOMER
#'
#' @return
#' @export
#'
#' @examples
#' res <- useHOMER("findMotifs.pl")
#' system(res)
#' 
useHOMER <- function(function.name){
  
  R_LIB=.libPaths()
  PATH1=Sys.getenv("PATH")
  
  homer_Lib=file.path(R_LIB,"ChipSeq/homer/bin")
  
  Sys.setenv(PATH=paste0(homer_Lib,":",PATH1)) 
  cmd1 <- Sys.which(function.name)[[1]]
  
  PATH2=Sys.getenv("PATH")
  
  cat(PATH1,"\n")
  cat(PATH2,"\n")
  
  return(cmd1)
}