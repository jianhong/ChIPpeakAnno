input.count.file <- "~/Aimin/DropboxUmass/Aimin/Project/umw_paul_kaufman/SmallRNA052219_133716585/all_count_se.txt"
output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/umw_paul_kaufman/SmallRNA052219_133716585"
fileName <- "all_count_se_with_gene_symbol"

count.table <- read.table(input.count.file,header = T)
rownames(count.table) <- count.table$Geneid

library(biomaRt)

mart= useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart = mart)

G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=count.table$Geneid,mart= mart)

getGeneSymbolOut <- function(count.table,G_list,output.file.dir,fileName) {
  count.data.1 <- data.frame(count.table,ensembl_gene_id = rownames(count.table))
  count.data.2 <- merge(count.data.1,G_list,by="ensembl_gene_id",sort = FALSE)
  count.data.2[which(count.data.2$hgnc_symbol == ""),]$hgnc_symbol <- paste0(count.data.2[which(count.data.2$hgnc_symbol == ""),]$ensembl_gene_id,"-No_match")
  write.table(count.data.2,file = file.path(output.file.dir,paste0(fileName,".txt")),append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F)
}

getGeneSymbolOut(count.table,G_list,output.file.dir,fileName)
