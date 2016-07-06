estLibSize <- function(bamfiles, index=bamfiles, ...){
    res <- mapply(function(f, i){
        if(suppressMessages(testPairedEndBam(f, index=i))){
            countBam(f, index=i, 
                     param=ScanBamParam(flag=scanBamFlag(isPaired = TRUE, 
                                                         isProperPair = TRUE,
                                                         isSecondaryAlignment = FALSE,
                                                         isNotPassingQualityControls = FALSE)))[, "records"]/2
        }else{
            countBam(f, index=i, 
                     param=ScanBamParam(flag=scanBamFlag(isPaired = FALSE, 
                                                         isUnmappedQuery = FALSE,
                                                         isSecondaryAlignment = FALSE,
                                                         isNotPassingQualityControls = FALSE)))[, "records"]
        }
    }, bamfiles, index, SIMPLIFY = TRUE)
    
    res
}