#!/usr/bin/env Rscript

input.file <- "/Users/aiminyan/Aimin/DropboxUmass/Aimin/Project/umw_castilla_lucio/CBFB RIP-seq.xlsx"

output.file.dir <- "/Users/aiminyan/Aimin/DropboxUmass/Aimin/Project/umw_castilla_lucio"

require(gdata)
df = read.xls (input.file, sheet = 1, header = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ChIPpeakAnno")

require(ChIPpeakAnno)

library(biomaRt)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","external_gene_name"),values=df$Ensemble.gene_id,mart= mart)

#count.data.2 <- merge(df,G_list,by="ensembl_gene_id")
#count.data.2[which(count.data.2$mgi_symbol %in% ""),]$mgi_symbol <- paste0(count.data.2[which(count.data.2$mgi_symbol %in% ""),]$ensembl_gene_id,"-No_match")

x_name <- "CBFB_RIP_seq_with_hgnc_symbol_external_gene_name.txt"
write.table(G_list,file = file.path(output.file.dir,x_name),append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T)

head(df)
tail(df)

save.image(file = file.path(output.file.dir,"CBFB_RIP_seq_with_hgnc_symbol_external_gene_name.RData"))
savehistory(file = file.path(output.file.dir,"CBFB_RIP_seq_with_hgnc_symbol_external_gene_name.Rhistory"))
