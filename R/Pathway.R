flatten_list = function(x)
{
  if (typeof(x) != "list")
  {
    return(x)
  }
  sapply(x, function(y) paste(y, collapse = " ; "))
}

writeTibble <- function(tibble.input, output.file.name = tempfile())
{
  if (!dir.exists(dirname(output.file.name)))
  {
    dir.create(dirname(output.file.name), recursive = TRUE)
  }
  tibble.input %>% mutate_all(funs(flatten_list)) %>% write.csv(output.file.name)
}

gmt.file <- "~/Aimin/DropboxUmass/Aimin/Project/Pathway_Julie/c7.all.v6.2.symbols.gmt.txt"

convertPathway <- function(input.pathway.file,output.pathway.file){

gmt.data <- gsa.read.gmt(gmt.file)
output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/Pathway_Julie"
XX <- tibble(genes=flatten_list(gmt.data$genesets),geneset=flatten_list(gmt.data$geneset.names),descriptions=flatten_list(gmt.data$geneset.descriptions))
writeTibble(XX, output.file.name = file.path(output.file.dir,paste0(x_name,"pahtway.csv")))

}

