test_addGeneIDs<-function(){
    data(annotatedPeak)
    x<-addGeneIDs(annotatedPeak[1:6,],"org.Hs.eg.db",c("symbol","omim"))
}