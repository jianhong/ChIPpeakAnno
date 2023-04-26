#' get the oligonucleotide frequency
#' 
#' Prepare the oligonucleotide frequency for given Markov order.
#' 
#' 
#' @param sequence The sequences packaged in DNAStringSet, DNAString object or
#' output of function \link{getAllPeakSequence}.
#' @param MarkovOrder Markov order.
#' @return A numeric vector.
#' @author Jianhong Ou
#' @seealso See Also as \code{\link{oligoSummary}}
#' @keywords misc
#' @importFrom Biostrings oligonucleotideFrequency
#' @export
#' @examples
#'     library(seqinr)
#'     library(Biostrings)
#'     oligoFrequency(DNAString("AATTCGACGTACAGATGACTAGACT"))
#' 
oligoFrequency <- function(sequence, MarkovOrder=3L){
    stopifnot(is.numeric(MarkovOrder))
    stopifnot(MarkovOrder>0)
    stopifnot("The seqinr package is required." = 
                  requireNamespace("seqinr", quietly = TRUE))
    if(is(sequence, "GRanges")){
        sequence <- sequence$sequence
    }
    if(is.character(sequence)){
      sequence <- DNAStringSet(sequence)
    }
    if(!inherits(sequence, c("DNAStringSet", "DNAString"))){
      stop("sequence must be an object of DNAStringSet or DNAString",
           "or output of getAllPeakSequence")
    }
    MarkovOrder <- as.integer(MarkovOrder)
    total <- sum(nchar(sequence))
    freqs <- lapply(unique(c(1, MarkovOrder, MarkovOrder+1)), function(m){
        of <- oligonucleotideFrequency(sequence,
                                       width = m,
                                       as.prob = FALSE)
        if(length(dim(of))==2) of <- colSums(of)
        of/total
    })
    unlist(freqs, recursive = FALSE)
}




#' Output a summary of consensus in the peaks
#' 
#' Calculate the z-scores of all combinations of oligonucleotide in a given
#' length by Markove chain.
#' 
#' 
#' @param sequence The sequences packaged in DNAStringSet, DNAString object or
#' output of function \link{getAllPeakSequence}.
#' @param oligoLength The length of oligonucleotide.
#' @param freqs Output of function \link{frequency}.
#' @param MarkovOrder The order of Markov chain.
#' @param quickMotif Generate the motif by z-score of not.
#' @param revcomp Consider both the given strand and the reverse complement
#' strand when searching for motifs in a complementable alphabet (ie DNA).
#' Default, FALSE.
#' @param maxsize Maximum allowed dataset size (in length of sequences).
#' @return A list is returned. \item{zscore}{A numeric vector. The z-scores of
#' each oligonucleotide.} \item{counts}{A numeric vector. The counts number of
#' each oligonucleotide.} \item{motifs}{a list of motif matrix.}
#' @author Jianhong Ou
#' @seealso See Also as \code{\link{frequency}}
#' @references van Helden, Jacques, Marcel li del Olmo, and Jose E.
#' Perez-Ortin. "Statistical analysis of yeast genomic downstream sequences
#' reveals putative polyadenylation signals." Nucleic Acids Research 28.4
#' (2000): 1000-1010.
#' @keywords misc
#' @export
#' @importFrom Biostrings DNAStringSet PDict vcountPDict DNAString 
#' pairwiseAlignment aligned
#' @importFrom utils adist combn
#' @importFrom stats hclust kmeans as.dendrogram nobs
#' @examples
#' 
#'     if(interactive() || Sys.getenv("USER")=="jianhongou"){
#'         data(annotatedPeak)
#'         library(BSgenome.Hsapiens.UCSC.hg19)
#'         library(seqinr)
#'         seq <- getAllPeakSequence(annotatedPeak[1:100], 
#'                      upstream=20, 
#'                      downstream=20, 
#'                      genome=Hsapiens)
#'         oligoSummary(seq)
#'     }
#' 
oligoSummary <- function(sequence, oligoLength=6L, 
                         freqs=NULL, MarkovOrder=3L, 
                         quickMotif=FALSE, revcomp=FALSE,
                         maxsize=100000){
    oligoLength <- as.integer(oligoLength)
    stopifnot(oligoLength>3&oligoLength<13)
    MarkovOrder <- as.integer(MarkovOrder)
    stopifnot(MarkovOrder>0&MarkovOrder<6)
    stopifnot(MarkovOrder < oligoLength-2)
    stopifnot("The seqinr package is required." = 
                  requireNamespace("seqinr", quietly = TRUE))
    if(is(sequence, "GRanges")){
        sequence <- sequence$sequence
    }else{
        if(is(sequence, "DNAStringSet")){
            sequence <- as.character(sequence)
        }else{
            if(!is.character(sequence)){
                stop("sequence must be an object of DNAStringSet or DNAString",
                     "or output of getAllPeakSequence")
            }
        }
    }
    
    sequence <- tolower(sequence)
    oligoWords <- seqinr::words(oligoLength)
    dict <- PDict(DNAStringSet(oligoWords))
    len <- length(sequence)
    sequence.tbl <- 1
    if(len > maxsize){
        sequence.tbl <- table(sequence)
        if(length(sequence.tbl)>maxsize){
            stop("oligoSummary can not handle such hugue dataset.")
        }
        cnt <- vcountPDict(dict, subject=DNAStringSet(names(sequence.tbl)), 
                           max.mismatch=0, min.mismatch=0, 
                           with.indels=FALSE, fixed=TRUE)
    }else{
        cnt <- vcountPDict(dict, subject=DNAStringSet(sequence), 
                           max.mismatch=0, min.mismatch=0, 
                           with.indels=FALSE, fixed=TRUE)
    }
    rownames(cnt) <- oligoWords
    cnt <- t(cnt)
    mergeRevcomp <- function(mat){
        coln <- colnames(mat)
        map <- c(a="t", c="g", g="c", t="a")
        revComp <- function(.ele){
            paste(map[rev(seqinr::s2c(.ele))], collapse="")
        }
        coln.rev <- sapply(coln, revComp)
        coln.id <- 1:length(coln)
        coln.rev.id <- match(coln.rev, coln)
        coln.ids <- apply(cbind(coln.id, coln.rev.id), 1, sort)
        coln.ids <- unique(t(coln.ids))
        mat1 <- mat[, coln.ids[, 1]] > 0
        mat2 <- mat[, coln.ids[, 2]] > 0
        mat1 <- colSums(mat1)
        mat2 <- colSums(mat2)
        coln <- ifelse(mat1>mat2, names(mat1), names(mat2))
        mat <- mat[, coln.ids[, 1]] + mat[, coln.ids[, 2]]
        colnames(mat) <- coln
        mat
    }
    if(revcomp) cnt <- mergeRevcomp(cnt)
    #cnt <- cnt>0
    mode(cnt) <- "logical"
    mode(cnt) <- "numeric"
    if(len > maxsize) cnt <- cnt * as.numeric(sequence.tbl)
    cntSum <- colSums(cnt)
    ## calculate Z score
    ## shuffle sequence based on 1 level markovmodel
    if(is.null(freqs)){
        freqs <- oligoFrequency(sequence, 
                                MarkovOrder=MarkovOrder)
    }
    namesFreqs <- unique(c("a", "c", "g", "t", 
                           seqinr::words(MarkovOrder),
                           seqinr::words(MarkovOrder+1)))
    names(freqs) <- tolower(names(freqs))
    stopifnot(all(namesFreqs %in% names(freqs)))
    f <- sapply(names(cntSum), function(.ele) {
        m1 <- substring(.ele, 1:(oligoLength-MarkovOrder),
                        (MarkovOrder+1):oligoLength)
        m0 <- substring(.ele, 2:(oligoLength-MarkovOrder),
                        (MarkovOrder+1):(oligoLength-1))
        prod(freqs[m1])/prod(freqs[m0])
    })
    if(len > maxsize){
        seqLen <- sapply(names(sequence.tbl), nchar, USE.NAMES=FALSE) - 
            oligoLength + 1
    }else{
        seqLen <- sapply(sequence, nchar, USE.NAMES=FALSE) - oligoLength + 1
    }
    seqLen[seqLen<0] <- 0
    f.m <- seqLen %*% t(f)
    f.m[f.m>1] <- 1
    if(len > maxsize){
        f.m <- f.m * as.numeric(sequence.tbl)
    }
    mu <- colSums(f.m)
    names(mu) <- names(cntSum)
    #     Kov <- sapply(names(cntSum), function(.ele){
    #         arr <- seqinr::s2c(.ele)
    #         k <- lapply(seq_len(length(arr)), function(.e){
    #             if(.e==1) return(TRUE)
    #             if(.e %% 2) 
    #                 return(all(arr[1:ceiling(.e/2)]==arr[ceiling(.e/2):.e]))
    #             return(all(arr[1:(.e/2)]==arr[(.e/2+1):.e]))
    #         })
    #         k <- as.numeric(k)
    #         j <- sapply(seq_len(length(arr)), function(.e){
    #             (1/freqs[arr[.e]])^.e
    #         })
    #         sum(k*j)
    #     })
    #std <- sqrt(mu * (2*Kov - 1 - (2*oligoLength - 1) * mu/length(sequence)))
    #std <- sqrt(length(sequence) * 
    #            (Kov/2-1+(2*oligoLength-1)/4^oligoLength) / 4^oligoLength)
    N <- len #length(sequence)
    std <- sqrt(mu/N * (1-mu/N) / N)
    zscore <- (cntSum - mu)/N/std
    #zscore <- (cntSum - mu)/std
    zscore[std==0] <- NA
    
    if(!quickMotif){
        return(list(zscore=zscore, counts=cntSum, expCnt=mu,
                    motifs=NA))
    }
    ## motif search
    seeds <- zscore[!is.na(zscore)]
    seeds <- sort(seeds, decreasing=TRUE)
    seeds <- seeds[1:min(50, length(seeds))]
    zscore.max <- max(zscore)
    seeds <- seeds[zscore.max - seeds <= zscore.max/2]
    zscore.cutoff <- c(NA, NA, 2.155, 2.66, 3.095, 3.49, 3.83, 4.1, 5, 5, 5, 5)
    seeds <- names(seeds)[seeds>zscore.cutoff[oligoLength]] #top50 only
    
    str2motif <- function(s){
        tmp <- matrix(c(1, 0, 0, 0,
                        0, 1, 0, 0,
                        0, 0, 1, 0,
                        0, 0, 0, 1), ncol=4, nrow=4, 
                      dimnames=list(c("A", "C", "G", "T"), 
                                    c("a", "c", "g", "t")))
        if(length(s)==1){
            return(tmp[, seqinr::s2c(s)])
        }
        ss <- table(s)
        if(length(ss)==1){
            return(tmp[, seqinr::s2c(names(ss))])
        }
        sss <- lapply(names(ss), function(.ele) 
            as.numeric(tmp[, seqinr::s2c(.ele)]))
        sss <- mapply(function(mat, times, coln) mat*times, 
                      sss, ss, SIMPLIFY = FALSE)
        consensusStr <- DNAString(names(ss)[1])
        for(i in 2:length(ss)){
            consensusStr <- aligned(
                    pairwiseAlignment(consensusStr, DNAString(names(ss)[i]), 
                                      type="local"),
                    degap=TRUE)
        }
        consensusStr <- tolower(as.character(consensusStr))
        idx <- regexpr(consensusStr, names(ss))
        tl <- max(idx)
        nc <- nchar(names(ss)[1])
        tmp0 <- rep(0, 4*(tl+nc-1))
        pcm <- mapply(function(.ele, .idx){
            left <- max(idx) - .idx
            this.tmp <- tmp0
            this.tmp[left*4+1:length(.ele)] <- .ele
            this.tmp
        }, sss, idx, SIMPLIFY=FALSE)
        pcm <- do.call(rbind, pcm)
        pcm <- colSums(pcm)
        pcm <- matrix(pcm, nrow=4, 
                      dimnames=list(c("A", "C", "G", "T")))
        pcm.colsums <- colSums(pcm)
        pcm.colsums.max <- max(pcm.colsums)
        pcm <- t(t(pcm)+(pcm.colsums.max - pcm.colsums)/4)
        pcm/colSums(pcm)
    }
    subgroupMotif <- function(.ele){
        .ele <- .ele[!is.na(.ele)]
        if(length(.ele)==0){
            return(NA)
        }
        if(length(.ele)==1){
            return(str2motif(.ele))
        }
        .cnt <- cnt[, .ele]
        .ord <- colSums(.cnt)
        .cnt <- .cnt[, order(-.ord)]
        if(length(.ele)==2){
            .ele <- colnames(.cnt)
            .ele <- c(rep(.ele[1], sum(.cnt[, 1])), 
                      rep(.ele[2], sum(.cnt[.cnt[, 2]>0 &
                                                (!(.cnt[, 1]>0 & 
                                                       .cnt[, 2]>0)), 2])))
            return(str2motif(.ele))
        }
        ## get max hits of any combinations
        ad <- adist(.ele, partial=TRUE)
        rownames(ad) <- colnames(ad) <- .ele
        all.comb <- sapply(1:min(length(.ele), 5), function(m){ 
            .comb <- combn(.ele, m, simplify=FALSE)
            .comb.ad <- sapply(.comb, function(.cmb){
                if(length(.cmb)>1){
                    .cmb <- combn(.cmb, 2)
                    all(ad[.cmb[1, ], .cmb[2, ]]<=2)
                }else{
                    TRUE
                }
            })
            .comb[.comb.ad]
        })
        all.comb <- unlist(all.comb, recursive=FALSE, use.names=FALSE)
        prob.mul.evt <- function(x){
            if(length(x)==1) return(x)
            p <- x[1] + x[2] - x[1]*x[2]
            Recall(c(p, x[-(1:2)]))
        }
        comb.f <- sapply(all.comb, function(.e) prob.mul.evt(f[.e]), 
                         simplify = TRUE,
                         USE.NAMES = FALSE)
        comb.f.m <- seqLen %*% t(comb.f)
        comb.f.m[comb.f.m>1] <- 1
        if(len > maxsize){
            comb.f.m <- comb.f.m * as.numeric(sequence.tbl)
        }
        comb.mu <- colSums(comb.f.m)
        rm(comb.f.m); gc(reset=TRUE)
        comb.std <- sqrt(comb.mu/N * (1-comb.mu/N)/N)
        comb.score <- sapply(all.comb, function(.e) sum(.cnt[, .e]),
                             simplify = TRUE,
                             USE.NAMES = FALSE)
        score <- (comb.score - comb.mu)/N/comb.std
        best.combs <- all.comb[score==max(score)]
        best.combs.score <- sapply(best.combs, function(.ele) sum(.ord[.ele]))
        best.comb <- best.combs[best.combs.score==max(best.combs.score)][[1]]
        .cnt <- .cnt[, colnames(.cnt) %in% best.comb, drop=FALSE]
        ## resort the column order
        .cnt <- .cnt[, order(.ord[colnames(.cnt)], decreasing=TRUE), drop=FALSE]
        .cnt.gt0 <- .cnt > 0
        .ele <- rep(colnames(.cnt)[1], sum(.cnt[, 1]))
        if(ncol(.cnt)>1){
            for(i in 2:ncol(.cnt)){
                .cnt.gt0[, i] <- .cnt.gt0[, i] & 
                    !(apply(.cnt.gt0[, 1:(i-1), drop=FALSE], 1, any) &
                          .cnt.gt0[, i])
                .ele <- c(.ele, rep(colnames(.cnt)[i], 
                                    sum(.cnt[.cnt.gt0[, i] , i])))
            }
        }
        str2motif(.ele)
    }
    if(length(seeds)>2){
        ### split the seeds by distance
        dist <- adist(seeds, partial=TRUE)
        hc <- hclust(dist(dist), method="average")
        len <- length(hc$height)
        if(len>2){
            hc.height.diff <- diff(hc$height)
            km <- kmeans(hc.height.diff, centers=2)
            idx <- km$cluster==km$cluster[length(km$cluster)]
            idx <- (length(idx):1)[which(rev(!idx))[1]]+1
            d <- cut(as.dendrogram(hc), hc$height[idx])
            lowers <- d$lower[sapply(d$lower, nobs)>=1]
            if(length(lowers)>0){
                subgroup <- lapply(lowers, function(.ele) 
                    seeds[as.numeric(labels(.ele))])
                subgroup <- subgroup[order(sapply(subgroup, function(.ele) 
                    max(zscore[.ele])), decreasing=TRUE)]
                motifs <- lapply(subgroup, subgroupMotif)
            }else{
                motifs <- NA ## this is impossible
            }
        }else{
            d <- as.dendrogram(hc)
            subgroup <- list()
            subgroup[[1]] <- seeds[as.numeric(labels(d[[1]]))]
            subgroup[[2]] <- seeds[as.numeric(labels(d[[2]]))]
            motifs <- lapply(subgroup, subgroupMotif)
        }
    }else{
        if(length(seeds)==2){
            consensus <- tolower(as.character(aligned(
                pairwiseAlignment(DNAString(seeds[1]), 
                                  DNAString(seeds[2]),
                                  type="local"),
                degap=TRUE)))
            if(nchar(consensus)>=oligoLength-2){
                motifs <- list(subgroupMotif(seeds))
            }else{
                motifs <- lapply(seeds[order(-zscore[seeds])], str2motif)
            }
        }else{
            motifs <- lapply(seeds[order(-zscore[seeds])], str2motif)
        }
    }
    
    return(list(zscore=zscore, counts=cntSum,
                motifs=motifs))
}
