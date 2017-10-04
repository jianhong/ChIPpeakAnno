xget <- function(x, envir, mode, ifnotfound=NA, inherits, 
                 output=c("all", "first")){
    output <- match.arg(output)
    y <- mget(x=x, envir=envir, mode=mode, 
              ifnotfound=ifnotfound, inherits=inherits)
    switch(output, 
           all=sapply(y, base::paste, collapse=";"),
           first=sapply(y, `[`, 1)
    )
}
