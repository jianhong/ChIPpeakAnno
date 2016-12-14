#' GetSampleInfo
#'
#' GetSampleInfo 
#' 
#' @param input.file 
#'  
#' 
#' @return
#' @export
#'
#' @examples
#' 
#' input.file="/Volumes/Bioinformatics$/2016/Danny/SampleID_INFO_ChIP.csv"
#'
#' re<-GetSampleInfo(input.file)
#' 
#' 
GetSampleInfo<- function(input.file) {
  
  re<-read.csv(input.file)
  re.c<-colnames(re)
  
  #print(re.c[2])
  #print(re.c[3])
  
  Ab.type<-unique(as.character(re[,colnames(re)==re.c[3]]))
  Cell.type<-unique(as.character(re[,colnames(re)==re.c[2]]))
  
  #print(Ab.type)
  #print(Cell.type)

  re1<-cbind(re[,1:3],paste0(re[,2],"_",re[,3]))
  
  colnames(re1)[4]="Cell_TF"
  
  Cell.Ab.type<-unique(as.character(re1[,4]))
  
  re2<-lapply(1:length(Cell.Ab.type),function(u,Cell.Ab.type,re1){
    
    x=Cell.Ab.type[u]
    z=re1
    
    ZZ<-as.character(z[which(z[,4]==x),1])
    ZZ
      
  },Cell.Ab.type,re1)

  names(re2)<-Cell.Ab.type
  
  #Get the combination of cell type and antibody with replicates
  
  re21<-re2[lapply(re2, length) == 1]
  
  re22<-re2[lapply(re2, length) == 2]
  
  re3<-list(re1=re1,re2=re2,re21=re21,re22=re22)
  
  return(re3)
  
}