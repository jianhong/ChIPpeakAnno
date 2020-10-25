#' Convert between the name of the organism annotation package ("OrgDb") and
#' the name of the organism.
#' 
#' Give a species name and return the organism annotation package name or give
#' an organism annotation package name then return the species name.
#' 
#' 
#' @param name The name of the organism annotation package or the species.
#' @return A object of character
#' @author Jianhong Ou
#' @keywords misc
#' @export
#' @importFrom utils adist
#' @examples
#' 
#'   egOrgMap("org.Hs.eg.db")
#'   egOrgMap("Mus musculus")
#' 
egOrgMap <- function(name){
    if(!is(name, "character"))
        stop("class of input organism should be character")
    organism <- c(
        "org.Ag.eg.db"="Anopheles gambiae",
        "org.At.eg.db"="Arabidopsis thaliana",
        "org.Bt.eg.db"="Bos taurus",
        "org.Ce.eg.db"="Caenorhabditis elegans",
        "org.Cf.eg.db"="Canis lupus familiaris",
        "org.Dm.eg.db"="Drosophila melanogaster",
        "org.Dr.eg.db"="Danio rerio",
        "org.EcK12.eg.db"="Escherichia coli str. K12a",
        "org.EcSakai.eg.db"="Escherichia coli O157:H7 str. Sakai",
        "org.Gg.eg.db"="Gallus gallus",
        "org.Hs.eg.db"="Homo sapiens",
        "org.Mm.eg.db"="Mus musculus",
        "org.Mmu.eg.db"="Macaca mulatta",
        "org.Pf.plasmo.db"="Plasmodium falciparum",
        "org.Pt.eg.db"="Pan troglodytes",
        "org.Rn.eg.db"="Rattus norvegicus",
        "org.Sc.sgd.db"="Saccharomyces cerevisiae",
        "org.Sco.eg.db"="Streptomyces coelicolor",
        "org.Ss.eg.db"="Sus scrofa",
        "org.Tgondii.eg.db"="Toxoplasma gondii",
        "org.Xl.eg.db"="Xenopus laevis")
    if(name %in% names(organism)){
        return(organism[name])
    }else{
        if(name %in% organism)
            return(names(organism)[organism==name])
    }
    ##claculate the string distance, need utils package
    org <- as.character(c(names(organism), organism))
    dis <- adist(name, org)
    org <- org[dis==min(dis)]
    stop(paste("You input \"", name, "\". Do you mean \"", org, "\"?", sep=""))
}
