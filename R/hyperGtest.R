hyperGtest <- function(alltermcount,thistermcount, 
                       totaltermInGenome, totaltermInPeakList){
    pvalue = as.numeric()
    thistermtotal = as.numeric()
    for (i in 1:length(thistermcount$GOterm))
    {
        m = 
            as.numeric(alltermcount$GOcount[alltermcount$GOterm==
                                                thistermcount$GOterm[i]])
        q = as.numeric(thistermcount$GOcount[i])
        n = as.numeric(totaltermInGenome)
        k = as.numeric(totaltermInPeakList)
        ####### k - number of total GO terms in the peak list
        ####### n - number of GO terms total
        ####### q - number of this GO terms in the peak list
        ####### m - number of this GO terms total
        pvalue[i] = phyper(q-1, m, n-m, k, lower.tail = FALSE, log.p = FALSE)
        thistermtotal[i] = m
    }
    list(thisterm=thistermcount$GOterm, 
         thistermcount=thistermcount$GOcount,
         thistermtotal=thistermtotal, 
         pvalue=pvalue, 
         totaltermInPeakList=k, 
         totaltermInGenome=n)
}

