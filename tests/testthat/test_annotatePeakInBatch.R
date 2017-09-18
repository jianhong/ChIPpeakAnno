test_that("annotatePeakInBatch works not correct", {
    load(system.file("extdata", "annotatedPeaks.2.10.0.rds", 
                     ## this file should be rewrote for fix the new format
                     package="ChIPpeakAnno"))
    checkRes <- function(step, s){
        s <- s[!is.na(s$feature)]
        if(length(s$feature)>0){
            s <- as.data.frame(s)
            step <- step[rownames(s), ]
            step <- as.data.frame(step)
            for(i in c("peak", "feature", "start_position", 
                       "end_position", "distancetoFeature", 
                       "shortestDistance")){
                #if(any(step[,i]!=s[,i])) stop(paste(i, "is not identical!"))
                expect_equal(step[, i], s[, i], info=paste(i, "is not identical!"))
            }
        }
    }
    ##check one row
    myPeak = GRanges("chr1", 
                     IRanges(start = c(17208381), 
                             end = c(17208381), 
                             names=c("Site1")), 
                     strand = c('+'))
    feature = GRanges("chr1", 
                      IRanges(start = c(17066767, 17180899), 
                              end = c(17267729, 17180971), 
                              names =c("Site1", "Site2")),
                      strand = c('+'))
    suppressWarnings(annotation.step1 <- 
                         annotatePeakInBatch(myPeak, AnnotationData=feature, 
                                             output="overlapping",
                                             multiple=FALSE))
    checkRes(annotation.step1,annotation.s1)
    ##check one
    suppressWarnings(annotation.step2 <- 
                         annotatePeakInBatch(myPeak, 
                                             AnnotationData=feature[1,], 
                                             output="overlapping",
                                             multiple=FALSE))
    checkRes(annotation.step2,annotation.s2)
    suppressWarnings(annotation.step3 <- 
                         annotatePeakInBatch(myPeak, 
                                             AnnotationData=feature[2,], 
                                             output="overlapping"))
    #checkRes(annotation.step3,annotation.s3)
    suppressWarnings(annotation.step4 <- 
                         annotatePeakInBatch(myPeak, 
                                             AnnotationData=feature[2,], 
                                             output="both"))
    checkRes(annotation.step4,annotation.s4)
    suppressWarnings(annotation.step5 <- 
                         annotatePeakInBatch(myPeak, 
                                             AnnotationData=feature[2,]))
    checkRes(annotation.step5,annotation.s5)
    ##check select first
    suppressWarnings(annotation.step6 <- 
                         annotatePeakInBatch(myPeak, 
                                             AnnotationData=feature, 
                                             output="overlapping",
                                             select="first"))
    checkRes(annotation.step6,annotation.s6)
    suppressWarnings(annotation.step7 <- 
                         annotatePeakInBatch(myPeak, AnnotationData=feature, 
                                             output="overlapping",
                                             select="last"))
    checkRes(annotation.step7,annotation.s7)
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
    myPeak = GRanges("1", IRanges(start=c(2, 7, 10), end=c(6, 8, 10), 
                                  names=letters[1:3]), 
                     strand="*")
    feature = GRanges("1", IRanges(start=c(1,3, 6, 8, 2, 3, 1, 3, 1), 
                                   end=c(2, 3, 7, 9, 6, 5, 5, 7, 7), 
                                   names=LETTERS[1:9]), 
                      strand="+")
    suppressWarnings(annotation.step8 <- 
                         annotatePeakInBatch(myPeak, AnnotationData=feature))
    checkRes(annotation.step8,annotation.s8)
    
    suppressWarnings(annotation.step9 <- 
                         annotatePeakInBatch(myPeak, AnnotationData=feature, 
                                             output="overlapping"))
    checkRes(annotation.step9,annotation.s9)
    
    suppressWarnings(annotation.step10 <- 
                         annotatePeakInBatch(myPeak, AnnotationData=feature, 
                                             output="both"))
    checkRes(annotation.step10,annotation.s10)
    
    #contain negative strand
    feature = GRanges("1", 
                      IRanges(start=c(1,3, 6, 8, 2, 3, 1, 3, 1), 
                              end=c(2, 3, 7, 9, 6, 5, 5, 7, 7), 
                              names=LETTERS[1:9]), 
                      strand=c("-", rep("+", 5), "-", "+", "+"))
    suppressWarnings(annotation.step11 <- 
                         annotatePeakInBatch(myPeak, AnnotationData=feature))
    checkRes(annotation.step11,annotation.s11)
    
    #consider strand
    myPeak = GRanges("1", 
                     IRanges(start=c(2, 7, 10), 
                             end=c(6, 8, 10), 
                             names=letters[1:3]), 
                     strand=c("-","+","+"))
    suppressWarnings(annotation.step12 <- 
                         annotatePeakInBatch(myPeak, AnnotationData=feature))
    checkRes(annotation.step12,annotation.s12)
    
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
    myPeak = GRanges("1", 
                     IRanges(start=c(2, 7, 10, 17), 
                             end=c(6, 8, 12, 17), 
                             names=letters[1:4]), 
                     strand="*")
    feature = GRanges("1", 
                      IRanges(start=c(1, 3, 6, 8, 2, 3, 1, 
                                      3, 1, 11, 13, 19, 18),
                              end=c(2, 3, 7, 9, 6, 5, 5, 
                                    7, 7, 13, 15, 20, 20),
                              names=LETTERS[1:13]), 
                      strand=c("-", "+", "+", "-", "+", "+", "-", 
                               "+", "+", "+", "-", "+", "+"))
    suppressWarnings(annotation <- 
                         annotatePeakInBatch(myPeak, 
                                             AnnotationData=feature, 
                                             maxgap=1,
                                             output="upstream&inside"))
    expect_equal(sum(annotation$peak=="a"& !is.na(annotation$feature)), 8)
    expect_equal(sum(annotation$peak=="b"& !is.na(annotation$feature)), 4)
    expect_equal(sum(annotation$peak=="c"& !is.na(annotation$feature)), 2)
    expect_equal(sum(annotation$peak=="d"& !is.na(annotation$feature)), 1)
    suppressWarnings(annotation <- 
                         annotatePeakInBatch(myPeak, 
                                             AnnotationData=feature, 
                                             maxgap=2,
                                             output="upstream&inside"))
    expect_equal(sum(annotation$peak=="a"& !is.na(annotation$feature)), 8)
    expect_equal(sum(annotation$peak=="b"& !is.na(annotation$feature)), 5)
    expect_equal(sum(annotation$peak=="c"& !is.na(annotation$feature)), 2)
    expect_equal(sum(annotation$peak=="d"& !is.na(annotation$feature)), 3)
    suppressWarnings(annotation <- 
                         annotatePeakInBatch(myPeak, 
                                             AnnotationData=feature, 
                                             maxgap=1,
                                             output="upstream"))
    expect_equal(sum(annotation$peak=="a"& !is.na(annotation$feature)), 6)
    expect_equal(sum(annotation$peak=="b"& !is.na(annotation$feature)), 0)
    expect_equal(sum(annotation$peak=="c"& !is.na(annotation$feature)), 2)
    expect_equal(sum(annotation$peak=="d"& !is.na(annotation$feature)), 1)
    suppressWarnings(annotation <- 
                         annotatePeakInBatch(myPeak, 
                                             AnnotationData=feature, 
                                             maxgap=2,
                                             output="upstream"))
    expect_equal(sum(annotation$peak=="a"& !is.na(annotation$feature)), 6)
    expect_equal(sum(annotation$peak=="b"& !is.na(annotation$feature)), 1)
    expect_equal(sum(annotation$peak=="c"& !is.na(annotation$feature)), 2)
    expect_equal(sum(annotation$peak=="d"& !is.na(annotation$feature)), 3)
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
    suppressWarnings(annotation <- 
                         annotatePeakInBatch(myPeak, 
                                             AnnotationData=feature, 
                                             maxgap=1,
                                             output="inside&downstream"))
    expect_equal(sum(annotation$peak=="a"& !is.na(annotation$feature)), 8)
    expect_equal(sum(annotation$peak=="b"& !is.na(annotation$feature)), 5)
    expect_equal(sum(annotation$peak=="c"& !is.na(annotation$feature)), 2)
    expect_equal(sum(annotation$peak=="d"& !is.na(annotation$feature)), 0)
    suppressWarnings(annotation <- 
                         annotatePeakInBatch(myPeak, 
                                             AnnotationData=feature, 
                                             maxgap=2,
                                             output="inside&downstream"))
    expect_equal(sum(annotation$peak=="a"& !is.na(annotation$feature)), 9)
    expect_equal(sum(annotation$peak=="b"& !is.na(annotation$feature)), 6)
    expect_equal(sum(annotation$peak=="c"& !is.na(annotation$feature)), 2)
    expect_equal(sum(annotation$peak=="d"& !is.na(annotation$feature)), 0)
    suppressWarnings(annotation <- 
                         annotatePeakInBatch(myPeak, 
                                             AnnotationData=feature, 
                                             maxgap=1,
                                             output="downstream"))
    expect_equal(sum(annotation$peak=="a"& !is.na(annotation$feature)), 2)
    expect_equal(sum(annotation$peak=="b"& !is.na(annotation$feature)), 5)
    expect_equal(sum(annotation$peak=="c"& !is.na(annotation$feature)), 1)
    expect_equal(sum(annotation$peak=="d"& !is.na(annotation$feature)), 0)
    suppressWarnings(annotation <- 
                         annotatePeakInBatch(myPeak, 
                                             AnnotationData=feature, 
                                             maxgap=2,
                                             output="downstream"))
    expect_equal(sum(annotation$peak=="a"& !is.na(annotation$feature)), 3)
    expect_equal(sum(annotation$peak=="b"& !is.na(annotation$feature)), 6)
    expect_equal(sum(annotation$peak=="c"& !is.na(annotation$feature)), 1)
    expect_equal(sum(annotation$peak=="d"& !is.na(annotation$feature)), 0)
    #check upstream2downstream
    #          aa      bb
    #    12345678901234567890
    #    AAAAAAAAAAAAAA
    #       BBBBBBBB
    #       CC    DD
    #
    myPeak <- GRanges("1", 
                      IRanges(start=c(7, 15), end=c(8, 16), names=c("a","b")), 
                      strand="*")
    feature <- GRanges("1", 
                       IRanges(start=c(1, 4, 4, 10), 
                               end=c(14, 11, 5, 11), 
                               names=c("A", "B", "C", "D")), 
                       strand="+")
    suppressWarnings(annotation <- 
                         annotatePeakInBatch(myPeak, 
                                             AnnotationData=feature, 
                                             maxgap=1,
                                             output="upstreamORdownstream"))
    suppressWarnings(annotation1 <- 
                         annotatePeakInBatch(myPeak, 
                                             AnnotationData=feature, 
                                             maxgap=2,
                                             output="upstreamORdownstream"))
    expect_equal(length(annotation), 2)
    expect_equal(length(annotation1), 3)
    
    #check distance is correct. real data
    testPeaks <- GRanges("chr6", 
                         IRanges(c(114643757, 114644908, 114706674, 114790153, 
                                   114859230, 114860358, 114860615, 114862060, 
                                   114888901, 114910896), width=1), 
                         strand="+")
    data(TSS.mouse.GRCm38)
    suppressWarnings(annotation.step13 <- 
                         annotatePeakInBatch(testPeaks , 
                                             AnnotationData = TSS.mouse.GRCm38,
                                             FeatureLocForDistance = "geneEnd",
                                             output="both"))
    ##   identical(annotation.step13,annotation.s13) Faild in 2.10.0
    data(TSS.human.GRCh37)
    suppressWarnings(peakList <- 
                         toGRanges(system.file("extdata", "peaks_hg19.bed", 
                                               package="ChIPpeakAnno"),
                                   header=F))
    suppressWarnings(annotation.step14 <- 
                         annotatePeakInBatch(peakList, 
                                             AnnotationData=TSS.human.GRCh37, 
                                             FeatureLocForDistance="TSS", 
                                             PeakLocForDistance="middle", 
                                             maxgap=5000, select="first", 
                                             output="overlapping"))
    #checkRes(annotation.step14,annotation.s14) , 
    #difference: old 34 ENSG00000246463, new 34 ENSG00000181894
    
#    data(myPeakList)
#    data(TSS.human.NCBI36)
#    annotation.step15 = annotatePeakInBatch(myPeakList, 
#                                           AnnotationData=TSS.human.NCBI36)
#    checkRes(annotation.step15,annotation.s15)
##old
# RangedData with 3 rows and 9 value columns across 24 spaces
# space                 ranges |           peak      strand         feature start_position end_position insideFeature distancetoFeature
# <factor>              <IRanges> |    <character> <character>     <character>      <numeric>    <numeric>   <character>         <numeric>
#     1_19_116634269 ENSG00000221040        1 [116634634, 116634734] | 1_19_116634269           + ENSG00000221040      116622751    116622965    downstream             11883
# 11_11_71974217 ENSG00000222241        3 [ 71974527,  71974627] | 11_11_71974217           + ENSG00000222241       71959018     71959304    downstream             15509
# 11_13_43336708 ENSG00000222330        3 [ 43337016,  43337116] | 11_13_43336708           + ENSG00000222330       43331447     43331749    downstream              5569
# shortestDistance fromOverlappingOrNearest
# <numeric>              <character>
#     1_19_116634269 ENSG00000221040            11669             NearestStart
# 11_11_71974217 ENSG00000222241            15223             NearestStart
# 11_13_43336708 ENSG00000222330             5267             NearestStart
##new    
#     RangedData with 3 rows and 10 value columns across 24 spaces
#     space                 ranges | strand           peak         feature start_position end_position feature_strand insideFeature distancetoFeature
#     <factor>              <IRanges> |  <Rle>    <character>     <character>      <integer>    <integer>    <character>      <factor>         <numeric>
#         1_19_116634269 ENSG00000209170     chr1 [116634634, 116634734] |      * 1_19_116634269 ENSG00000209170      116622751    116622966              +    downstream             11883
#     11_13_43336708 ENSG00000210484     chr3 [ 43337016,  43337116] |      * 11_13_43336708 ENSG00000210484       43331447     43331751              +    downstream              5569
#     11_11_71974217 ENSG00000210995     chr3 [ 71974527,  71974627] |      * 11_11_71974217 ENSG00000210995       71959018     71959311              +    downstream             15509
#     shortestDistance fromOverlappingOrNearest
#     <integer>              <character>
#         1_19_116634269 ENSG00000209170            11668          NearestLocation
#     11_13_43336708 ENSG00000210484             5265          NearestLocation
#     11_11_71974217 ENSG00000210995            15216          NearestLocation
    
    peak <- GRanges(c("Chr01","Chr01","Chr01"),
                    IRanges(start=c(1650,2806,8361), 
                            end=c(1860,3006,8591),
                            names=c("peak1","peak2","peak3")),
                    strand=as.integer(1))

    GENOME <- GRanges("Chr01", 
                      IRanges(start=c(1660, 2906, 8391),
                              end=c(2502, 6646, 8860),
                              names=c("Potri.001G000100",
                                      "Potri.001G000200",
                                      "Potri.001G000300")),
                      strand=c("-", "-", "+"))
    suppressWarnings(annotation.step16 <- 
                         annotatePeakInBatch(peak, 
                                             AnnotationData=GENOME, 
                                             PeakLocForDistance ="middle",
                                             output="both"))
##    identical(annotation.step16,annotation.s16)
## there is a bug for 2.10.0 when calculating the distancetoFeature for negative strand.
    
    if(interactive()){
        PeakLocForDistance <- c("start","middle","end")
        FeatureLocForDistance <- c("TSS","middle","start","end", "geneEnd")
        select <- c("all", "first", "last", "arbitrary")
        output <- c("nearestLocation", "overlapping", "both", 
                    "shortestDistance", "upstream&inside", 
                    "inside&downstream", "upstream", "downstream", 
                    "upstreamORdownstream", "inside") 
        ##, "shortestDistance" is not in version 2.10.0
        comb <- do.call(expand.grid, 
                        list(PeakLocForDistance, 
                             FeatureLocForDistance, 
                             select, 
                             output))
        annotation.step17 <- apply(comb, 1, function(.ele) 
            annotatePeakInBatch(peak, AnnotationData=GENOME, 
                                PeakLocForDistance=.ele[1],
                                FeatureLocForDistance=.ele[2],
                                select=.ele[3],
                                output=.ele[4])) 
        ## CPU time long. maybe we can omit select option. 
    }
    
#annotation.s17 <- lapply(annotation.s17, function(.ele) {.ele$distancetoFeature <- abs(.ele$distancetoFeature); .ele})
#annotation.step17 <- lapply(annotation.step17, function(.ele) {.ele$distancetoFeature <- abs(.ele$distancetoFeature); .ele})
#checkRes(annotation.step17,annotation.s17)
    
    #check annotation by txdb
    txdb_file <- system.file("extdata", "Biomart_Ensembl_sample.sqlite",
                             package="GenomicFeatures")
    TxDb <- loadDb(txdb_file)
    exons <- unique(exons(TxDb, columns="exon_name"))
    names(exons) <- exons$exon_name
    fiveUTRs <- unique(unlist(fiveUTRsByTranscript(TxDb)))
    names(fiveUTRs) <- fiveUTRs$exon_name
    threeUTRs <- unique(unlist(threeUTRsByTranscript(TxDb)))
    names(threeUTRs) <- threeUTRs$exon_name
    gaps <- gaps(exons)
    gaps <- gaps[strand(gaps)!="*"]
    gaps.anno <- annotatePeakInBatch(gaps, 
                                     AnnotationData=exons, 
                                     output="inside")
    expect_true(all(is.na(gaps.anno$feature)))
    start(gaps) <- start(gaps)+2
    end(gaps) <- end(gaps)-2
    gaps.anno <- annotatePeakInBatch(exons, 
                                     AnnotationData=gaps, 
                                     output="overlapping", 
                                     maxgap=0,
                                     ignore.strand = FALSE)
    gaps.anno.1 <- gaps.anno[!is.na(gaps.anno$feature)]
    expect_equal(length(gaps.anno.1), 0)
    gaps.anno <- annotatePeakInBatch(exons, 
                                     AnnotationData=gaps, 
                                     output="overlapping", 
                                     maxgap=3,
                                     ignore.strand = FALSE)
    gaps.anno.2 <- gaps.anno[!is.na(gaps.anno$feature)]
    expect_equal(length(gaps.anno.2), 2*length(exons))
    exons.anno <- annotatePeakInBatch(exons,
                                      AnnotationData=exons, 
                                      output="inside")
    expect_equal(exons.anno$peak, exons.anno$feature)
    exons.anno <- annotatePeakInBatch(exons, 
                                      AnnotationData=exons, 
                                      output="nearestLocation")
    expect_equal(exons.anno$peak, exons.anno$feature)
    threeUTRs.anno <- annotatePeakInBatch(threeUTRs, 
                                          AnnotationData=exons, 
                                          output="downstream", maxgap=1000)
    fiveUTRs.anno <- annotatePeakInBatch(fiveUTRs, 
                                         AnnotationData=exons, 
                                         output="upstream", maxgap=1000)
    threeUTRs.anno <- annotatePeakInBatch(threeUTRs, 
                                          AnnotationData=exons, 
                                          output="inside&downstream")
    fiveUTRs.anno <- annotatePeakInBatch(fiveUTRs, 
                                         AnnotationData=exons, 
                                         output="upstream&inside")
    expect_equal(threeUTRs.anno$peak, threeUTRs.anno$feature)
    expect_equal(fiveUTRs.anno$peak, fiveUTRs.anno$feature)
})
