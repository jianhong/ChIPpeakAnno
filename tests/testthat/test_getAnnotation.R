test_that("getAnnotation works not correct", {
    #GRanges
    if(interactive()){
        mart<-useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")
        Annotation = getAnnotation(mart, featureType="TSS")
    }
})