
# Filtering list based on certain pattern of the names of list elements

FL <- function(re.bed, patt) {
  
  re.peaks.only.bed<-lapply(1:length(re.bed),function(u,re.bed,patt){
    
    tmp=re.bed
    x=tmp[[u]]
    path_name=dirname(x)
    file_name=basename(x)
    
    #pos_summits=regexpr('summits',file_name)
    pos_peaks=regexpr(patt,file_name)
    
    if(pos_peaks>0){
      y<-x
    }else{
      y<-NULL  
    }
    
    y  
    
  },re.bed,patt)
  
  re.peaks.only.bed.2<-re.peaks.only.bed[lapply(re.peaks.only.bed,length) > 0]
  
  names(re.peaks.only.bed.2)=unlist(lapply(1:length(re.peaks.only.bed.2),function(u,re.peaks.only.bed.2){
    tmp=re.peaks.only.bed.2
    x=tmp[[u]]
    path_name=dirname(x)
    file_name=basename(x)
    file_name
  },re.peaks.only.bed.2))
  
  return(re.peaks.only.bed.2)
}

list_to_df <- function(list_for_df) {
  list_for_df <- as.list(list_for_df)
  
  nm <- names(list_for_df)
  
  if (is.null(nm)) 
    nm <- seq_along(list_for_df)
  
  ID<-sapply(1:length(nm),function(u,nm){
    x=nm[u]
    pos=regexpr("\\.",x)
    pos1=pos-1
    y<-substr(x,1,pos1)
    y
  },nm)
  
  df <- data.frame(ID=ID,name = nm, stringsAsFactors = FALSE)
  df$value <- unname(list_for_df)
  df
}
