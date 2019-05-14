# function to get the overlapping region 
intersectBED <- function(inputGRangeList, byBase = TRUE, output = c("BED", "WIG"), outputFileNamePrefix = "common"){
  
  
  
}

uniqueBED <- function(inputGRangeList, byBase = TRUE, output  = c("BED", "WIG"), outputFileNamePrefix = c("file1Unique", "file2Unique"){
  
  
}

# The updated function for drawing heatmap around TSS that you mentioned to me previously

# Plot Jarcard heatmap with flexible parameters: 
plotHeatMapOfJI <- function(jarcardIndexMatrix, …){
  
  
}

plotAggregatedSignalsAroundPeaks <- function(peakGRangeList, upstream = 5000, downstream = 3000, combineSamples = TRUE, combineSamplesBy = c(mean, median), …)
       
# For example, NADs from NADfinder can be used as domainGRangesWithSRs,  NADs from MACs, normR and EDD can be used as domainGRangelistWithoutSRs, the output should be NSRs (NAD splitting regions) in this example. We name DSRs to be more general.                     
getDomainSplittingRegions <- function(domainGRangesWithSRs, domainGRangelistWithoutSRs)

