ratio.zScore <- function(A, B, background=1){
    stopifnot(length(A)==length(B))
    r <- log2(A+background) - log2(B+background)
    ## check distribution?
    ## normal distribution
    pop_sd <- sd(r) * seqrt((length(r)-1)/length(r))
    pop_mean <- mean(r)
    z <- (r - pop_mean)/pop_sd
    z
}