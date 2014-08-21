test_annotatePeakInBatch<-function(){
    ##check one row
    myPeak = RangedData(IRanges(start = c(17208381), end = c(17208381), names=c("Site1")),space = c("chr1"),strand = c('+'))
    feature = RangedData(IRanges(start = c(17066767, 17180899), end = c(17267729, 17180971), names =c("Site1", "Site2")),space = c("chr1"),strand = c('+'))
    suppressWarnings(annotation.step1 <- annotatePeakInBatch(myPeak, AnnotationData=feature, output="overlapping", maxgap=0, multiple=FALSE))
    ##check one
    suppressWarnings(annotation.step2 <- annotatePeakInBatch(myPeak, AnnotationData=feature[1,], output="overlapping", maxgap=0, multiple=FALSE))
    suppressWarnings(annotation.step3 <- annotatePeakInBatch(myPeak, AnnotationData=feature[2,], output="overlapping"))
    suppressWarnings(annotation.step4 <- annotatePeakInBatch(myPeak, AnnotationData=feature[2,], output="both"))
    suppressWarnings(annotation.step5 <- annotatePeakInBatch(myPeak, AnnotationData=feature[2,]))
    ##check select first
    suppressWarnings(annotation.step6 <- annotatePeakInBatch(myPeak, AnnotationData=feature, output="overlapping", maxgap=0, select="first"))
    suppressWarnings(annotation.step7 <- annotatePeakInBatch(myPeak, AnnotationData=feature, output="overlapping", maxgap=0, select="last"))
    ## algorithm
    ##
    #       aaaaabb c
    #      1234567890
    #      AABB CCDD 
    #       EEEEE
    #        FFF
    #      GGGGG
    #        HHHHH
    #      IIIIIII
    #ignore strand
    myPeak = RangedData(IRanges(start=c(2, 7, 10), end=c(6, 8, 10), names=letters[1:3]), space="1", strand="*")
    feature = RangedData(IRanges(start=c(1,3, 6, 8, 2, 3, 1, 3, 1), end=c(2, 3, 7, 9, 6, 5, 5, 7, 7), names=LETTERS[1:9]), space="1", strand="+")
    suppressWarnings(annotation.step8 <- annotatePeakInBatch(myPeak, AnnotationData=feature))
    
    suppressWarnings(annotation.step9 <- annotatePeakInBatch(myPeak, AnnotationData=feature, output="overlapping"))
    
    suppressWarnings(annotation.step10 <- annotatePeakInBatch(myPeak, AnnotationData=feature, output="both"))
    
    #contain negative strand
    feature = RangedData(IRanges(start=c(1,3, 6, 8, 2, 3, 1, 3, 1), end=c(2, 3, 7, 9, 6, 5, 5, 7, 7), names=LETTERS[1:9]), space="1", strand=c("-", rep("+", 5), "-", "+", "+"))
    suppressWarnings(annotation.step11 <- annotatePeakInBatch(myPeak, AnnotationData=feature))
    
    #consider strand
    myPeak = RangedData(IRanges(start=c(2, 7, 10), end=c(6, 8, 10), names=letters[1:3]), space="1", strand=c("-","+","+"))
    suppressWarnings(annotation.step12 <- annotatePeakInBatch(myPeak, AnnotationData=feature))
    
    #check upsteam&inside and inside&downstream
    #       aaaaabb ccc    d
    #      12345678901234567890
    #     -AA    -DD
    #       +BB+CC  +JJJ
    #      +EEEEE     -KKK
    #       +FFF            +LL
    #     -GGGGG           +MMM
    #       +HHHHH
    #     +IIIIIII
    
    #check distance is correct. real data
    testPeaks <- RangedData(IRanges(c(114643757, 114644908, 114706674, 114790153, 114859230, 114860358, 114860615, 114862060, 114888901, 114910896), width=1), space="chr6", strand="+")
    data(TSS.mouse.GRCm38)
    suppressWarnings(annotation.step13 <- annotatePeakInBatch(testPeaks , AnnotationData = TSS.mouse.GRCm38, FeatureLocForDistance = "geneEnd", output="both"))
    
    data(TSS.human.GRCh37)
    suppressWarnings(peakList <- BED2RangedData(system.file("extdata", "peaks_hg19.bed", 
                                                            package="ChIPpeakAnno"),
                                                header=F))
    suppressWarnings(annotation.step14 <- annotatePeakInBatch(peakList, AnnotationData=TSS.human.GRCh37, FeatureLocForDistance="TSS", PeakLocForDistance="middle", maxgap=5000, select="first", output="overlapping"))
    
    data(myPeakList)
    data(TSS.human.NCBI36)
    annotation.step15 = annotatePeakInBatch(myPeakList, AnnotationData=TSS.human.NCBI36)
    
    
    peak <- RangedData(IRanges(start=as.numeric(as.character(c(1650,2806,8361))), end=as.numeric(as.character(c(1860,3006,8591)))), space=c("Chr01","Chr01","Chr01"),strand=as.integer(1))
    row.names(peak)<-c("peak1","peak2","peak3")
    genome <- as.data.frame(rbind(c("Chr01","phytozome8_0",      "gene",1660,2502, ".",     "-",     ".","Potri.001G000100"),c("Chr01", "phytozome8_0",     "gene",2906,6646,    ".",        "-",     ".","Potri.001G000200"),c("Chr01", "phytozome8_0",   "gene",8391,8860, ".",     "+",     ".","Potri.001G000300")))
    suppressWarnings(GENOME<-GFF2RangedData(genome))
    row.names(GENOME)<-genome[,9]
#    suppressWarnings(annotation.step16 <- annotatePeakInBatch(peak, AnnotationData=GENOME, PeakLocForDistance ="middle",output="both"))
    
    PeakLocForDistance <- c("start","middle","end")
    FeatureLocForDistance <- c("TSS","middle","start","end", "geneEnd")
    select <- c("all", "first", "last", "arbitrary")
    output <- c("nearestStart", "overlapping", "both") ##, "shortestDistance" is not in version 2.10.0
    comb <- do.call(expand.grid, list(PeakLocForDistance, FeatureLocForDistance, select, output))
    annotation.step17 <- apply(comb, 1, function(.ele) 
        annotatePeakInBatch(peak, AnnotationData=GENOME, 
                            PeakLocForDistance=.ele[1],
                            FeatureLocForDistance=.ele[2],
                            select=.ele[3],
                            output=.ele[4])) ## CPU time long. maybe we can omit select option. 
    
    load(system.file("extdata", "annotatedPeaks.2.10.0.rds", 
                     package="ChIPpeakAnno"))
checkIdentical(annotation.step1,annotation.s1)
checkIdentical(annotation.step2,annotation.s2)
checkIdentical(annotation.step3,annotation.s3)
checkIdentical(annotation.step4,annotation.s4)
checkIdentical(annotation.step5,annotation.s5)
checkIdentical(annotation.step6,annotation.s6)
checkIdentical(annotation.step7,annotation.s7)
checkIdentical(annotation.step8,annotation.s8)
checkIdentical(annotation.step9,annotation.s9)
checkIdentical(annotation.step10,annotation.s10)
checkIdentical(annotation.step11,annotation.s11)
checkIdentical(annotation.step12,annotation.s12)
 ##   identical(annotation.step13,annotation.s13) Faild in 2.10.0
checkIdentical(annotation.step14,annotation.s14) 
    annotation.step15 <- annotation.step15[match(annotation.s15$peak, annotation.step15$peak),]
checkIdentical(annotation.step15,annotation.s15)
##    identical(annotation.step16,annotation.s16)
## there is a bug for 2.10.0 when calculating the distancetoFeature for negative strand.
annotation.s17 <- lapply(annotation.s17, function(.ele) {.ele$distancetoFeature <- abs(.ele$distancetoFeature); .ele})
annotation.step17 <- lapply(annotation.step17, function(.ele) {.ele$distancetoFeature <- abs(.ele$distancetoFeature); .ele})
checkIdentical(annotation.step17,annotation.s17)
}