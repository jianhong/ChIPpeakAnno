test_getAnnotation<-function(){
    #GRanges
    mart<-useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")
    Annotation = getAnnotation(mart, featureType="TSS", output="GRanges")
}