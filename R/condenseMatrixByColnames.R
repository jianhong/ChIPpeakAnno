##condense matrix by colnames
# parameters
# @mx    	a matrix to be condensed
# @iname	index colname
# @sep		separator
# @cnt		add count columns?
# output	index colname, condensed colnames with count columns.
# Usage:
#	a<-matrix(c(rep(rep(1:5,2),2),rep(1:10,2)),ncol=4)
#	colnames(a)<-c("con.1","con.2","index.1","index.2")
#	condenseMatrixByColnames(a,"con.1")
#	condenseMatrixByColnames(a,2)

condenseMatrixByColnames <- function(mx,iname,sep=";",cnt=FALSE){
    if(!is(mx, "matrix")) stop("mx must be matrix\n")
    if(length(iname)!=1) stop("iname must be single colname\n")
    m_cname<-colnames(mx)
    if(is(iname, "numeric") &iname<=length(m_cname)) iname<-m_cname[iname]
    cnames<-m_cname[which(m_cname!=iname)]
    if(length(m_cname)==length(cnames)) 
        stop("the colum name specified for condense does not exist")
    m_split<-split(mx[,cnames],mx[,iname])
    colN<-length(cnames)
    m_list<-lapply(m_split,function(.ele){
        x<-apply(matrix(.ele,nrow=colN,byrow=TRUE),1,base::paste,collapse=sep)
        if(cnt) {
            unlist(lapply(x, function(w){
                tmp<-unique(as.character(unlist(strsplit(w,sep))))
                c(paste(tmp,collapse=sep),length(tmp))
            }))
        }else{
            unlist(lapply(x, function(w){
                tmp<-unique(as.character(unlist(strsplit(w,sep))))
                paste(tmp,collapse=sep)}))
        }
    })
    m_dat<-as.data.frame(do.call(rbind,m_list))
    m_dat$index<-rownames(m_dat)
    cnames.cnt<-c()
    for(i in cnames){
        if(cnt) cnames.cnt<-c(cnames.cnt,i,paste(i,"count",sep="."))
        else cnames.cnt<-c(cnames.cnt,i)
    }
    colnames(m_dat)<-c(cnames.cnt,iname)
    m_dat<-m_dat[,c(iname,cnames.cnt)]
    m_dat
}
