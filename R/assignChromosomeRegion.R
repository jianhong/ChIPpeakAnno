assignChromosomeRegion <-
function(peaks.RD, exon, TSS, utr5, utr3, proximal.promoter.cutoff=1000, immediate.downstream.cutoff=1000)
{
	ann.peaks = annotatePeakInBatch(peaks.RD, AnnotationData = TSS)
	ann.peaks = ann.peaks[!is.na(ann.peaks$distancetoFeature),]
	upstream = ann.peaks[ ann.peaks$insideFeature=="upstream" | (ann.peaks$distancetoFeature<0 & ann.peaks$insideFeature =="overlapStart" & abs(ann.peaks$distancetoFeature)>ann.peaks$shortestDistance ) | ann.peaks$insideFeature=="includeFeature" | (ann.peaks$distancetoFeature>=0 & ann.peaks$insideFeature =="overlapStart" & ann.peaks$distancetoFeature ==ann.peaks$shortestDistance),]
 
	proximal.promoter.n = dim(upstream[upstream$distancetoFeature>= - proximal.promoter.cutoff | upstream$shortestDistance<= proximal.promoter.cutoff, ])[1]
	enhancer.n = dim(upstream)[1] - proximal.promoter.n

	downstream = ann.peaks[ann.peaks$insideFeature =="downstream",]
	immediateDownstream.n = dim(downstream[downstream$distancetoFeature <= immediate.downstream.cutoff,])[1]
	enhancer.n = enhancer.n +  dim(downstream[downstream$distancetoFeature > immediate.downstream.cutoff,])[1]

	 inside.peaks =  ann.peaks[ann.peaks$insideFeature =="inside" | ann.peaks$insideFeature =="overlapEnd" |  (ann.peaks$insideFeature =="overlapStart" & ann.peaks$distancetoFeature >=0 & ann.peaks$distancetoFeature != ann.peaks$shortestDistance) | (ann.peaks$insideFeature =="overlapStart" & ann.peaks$distancetoFeature <0 & abs(ann.peaks$distancetoFeature) ==ann.peaks$shortestDistance) ,]

	ann.utr5.peaks = annotatePeakInBatch(inside.peaks, AnnotationData = utr5)

	proximal.promoter.n = proximal.promoter.n + dim(ann.utr5.peaks[ann.utr5.peaks$insideFeature =="upstream",])[1]

	utr5.n = dim(ann.utr5.peaks[ann.utr5.peaks$insideFeature %in% c("includeFeature" , "inside") | (ann.utr5.peaks$insideFeature =="overlapStart"  & ann.utr5.peaks$distancetoFeature >=0 & ann.utr5.peaks$distancetoFeature != ann.utr5.peaks$shortestDistance)  | (ann.utr5.peaks$insideFeature =="overlapStart" & ann.utr5.peaks$distancetoFeature <0 & abs(ann.utr5.peaks$distancetoFeature) ==ann.utr5.peaks$shortestDistance)  | (ann.utr5.peaks$insideFeature =="overlapEnd" & ann.utr5.peaks$strand=="+" & abs(start(ann.utr5.peaks)-ann.utr5.peaks$end_position) >= (end(ann.utr5.peaks)-ann.utr5.peaks$end_position)) | (ann.utr5.peaks$insideFeature =="overlapEnd" & ann.utr5.peaks$strand=="-" & abs(end(ann.utr5.peaks)-ann.utr5.peaks$start_position) >= abs(start(ann.utr5.peaks)-ann.utr5.peaks$start_position )), ] )[1]

	proximal.promoter.n = proximal.promoter.n +  dim(ann.utr5.peaks[(ann.utr5.peaks$insideFeature =="overlapStart"  & ann.utr5.peaks$distancetoFeature >=0 & ann.utr5.peaks$distancetoFeature == ann.utr5.peaks$shortestDistance)  | (ann.utr5.peaks$insideFeature =="overlapStart" & ann.utr5.peaks$distancetoFeature <0 & abs(ann.utr5.peaks$distancetoFeature) !=ann.utr5.peaks$shortestDistance),])[1]

	downstream.utr5 = ann.utr5.peaks[ann.utr5.peaks$insideFeature =="downstream" | (ann.utr5.peaks$insideFeature =="overlapEnd" & ann.utr5.peaks$strand=="+" & abs(start(ann.utr5.peaks)-ann.utr5.peaks$end_position) < (end(ann.utr5.peaks)-ann.utr5.peaks$end_position)) | (ann.utr5.peaks$insideFeature =="overlapEnd" & ann.utr5.peaks$strand=="-" & abs(end(ann.utr5.peaks)-ann.utr5.peaks$start_position) < abs(start(ann.utr5.peaks)-ann.utr5.peaks$start_position )), ] 

	ann.utr3.peaks = annotatePeakInBatch(downstream.utr5, AnnotationData = utr3)

	utr3.n = dim(ann.utr3.peaks[ann.utr3.peaks$insideFeature %in% c("includeFeature" , "overlapStart", "overlapEnd", "inside"),])[1]

	rest.peaks = ann.utr3.peaks[ann.utr3.peaks$insideFeature %in% c("downstream", "upstream"),]

	ann.rest.peaks = annotatePeakInBatch(rest.peaks, AnnotationData = exon)

	intron.n = dim(ann.rest.peaks[ann.rest.peaks$insideFeature %in% c("downstream", "upstream"),])[1]
	exon.n =dim(ann.rest.peaks)[1] - intron.n

	total = dim(peaks.RD)[1]/100

	 list( "Exon" =exon.n/total, "Intron"=intron.n/total, "5UTR" = utr5.n/total, "3UTR" = utr3.n/total, "Proximal Promoter%"= proximal.promoter.n/total, "Immediate Downstream" = immediateDownstream.n/total, "Enhancer" = enhancer.n/total)
}
