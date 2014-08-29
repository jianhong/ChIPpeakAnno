test_peaksNearBDP<-function(){
    peaks <- BED2RangedData(system.file("extdata", "test.bed", 
                                           package="ChIPpeakAnno"),
                               header=F)
    annotation = GFF2RangedData(system.file("extdata", "ref.gtf", 
                                            package="ChIPpeakAnno"),
                                header=F, comment.char="#", sep="\t")
    peaksNearBDP(peaks, AnnotationData = annotation, MaxDistance = 3000, 
                 PeakLocForDistance = "middle", FeatureLocForDistance = "TSS")
}