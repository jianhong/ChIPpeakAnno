#' estimate the library size
#' 
#' estimate the library size of bam files
#' 
#' 
#' @param bamfiles The file names of the 'BAM' ('SAM' for asBam) files to be
#' processed.
#' @param index The names of the index file of the 'BAM' file being processed;
#' this is given without the '.bai' extension.
#' @param \dots Not used.
#' @return numberic vector
#' @author Jianhong Ou
#' @keywords misc
#' @export
#' @importFrom Rsamtools testPairedEndBam countBam ScanBamParam scanBamFlag
#' @examples
#' 
#'     if(interactive() || Sys.getenv("USER")=="jianhongou"){
#'         path <- system.file("extdata", "reads", package="MMDiffBamSubset")
#'         if(file.exists(path)){
#'             WT.AB2 <- file.path(path, "WT_2.bam")
#'             Null.AB2 <- file.path(path, "Null_2.bam")
#'             Resc.AB2 <- file.path(path, "Resc_2.bam")
#'             estLibSize(c(WT.AB2, Null.AB2, Resc.AB2))
#'         }
#'     }
#' 
estLibSize <- function(bamfiles, index=bamfiles, ...){
    res <- mapply(function(f, i){
        if(suppressMessages(testPairedEndBam(f, index=i))){
            countBam(f, index=i, 
                     param=ScanBamParam(
                         flag=scanBamFlag(
                             isPaired = TRUE, 
                             isProperPair = TRUE,
                             isSecondaryAlignment = FALSE,
                             isNotPassingQualityControls = FALSE)))[, 
                                                                    "records"]/2
        }else{
            countBam(f, index=i, 
                     param=ScanBamParam(
                         flag=scanBamFlag(isPaired = FALSE, 
                                          isUnmappedQuery = FALSE,
                                          isSecondaryAlignment = FALSE,
                                          isNotPassingQualityControls = 
                                              FALSE)))[, "records"]
        }
    }, bamfiles, index, SIMPLIFY = TRUE)
    
    res
}
