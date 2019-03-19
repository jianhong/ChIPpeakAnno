#!/usr/bin/env Rscript

output.file.dir <- "~/Aimin/DropboxUmass/Aimin/Project/Simone_Tessa"

input.file <- "~/Aimin/DropboxUmass/Aimin/Project/Simone_Tessa/GSE78220_PatientFPKM.xlsx"

rnaSeq.data <- as.data.frame(readxl::read_xlsx(input.file))
colnames(rnaSeq.data) <- tools::file_path_sans_ext(colnames(rnaSeq.data))

Non_responding <- paste0("Pt",c(1,7,10,11,12,14,16,17,20,22,23,25,29,30,31,32,36))
Responding <- paste0("Pt",c(2,3,4,5,6,8,9,13,15,18,19,21,24,26,"27A","27B",28,33,34,35,37,38))

length(Non_responding)
length(Responding)

gse78220 <- getGEO(filename='~/Aimin/DropboxUmass/Aimin/Project/Simone_Tessa/GSE78220_series_matrix.txt.gz')
phenoData(gse78220)

index <- match(Non_responding,colnames(rnaSeq.data))
index1 <- index[which(!is.na(index))]

rnaSeq.data.non.responding <- rnaSeq.data[,index1]
rnaSeq.data.responding <- rnaSeq.data[,-c(1,index1)]

rownames(rnaSeq.data.non.responding) <-  rnaSeq.data$Gene
rownames(rnaSeq.data.responding) <-  rnaSeq.data$Gene

data.non <- cbind(colnames(rnaSeq.data.non.responding),rep("non_responding",dim(rnaSeq.data.non.responding)[2]))
data.re <- cbind(colnames(rnaSeq.data.responding),rep("responding",dim(rnaSeq.data.responding)[2]))

cel.file.sample.infor.no.2<-as.data.frame(rbind(data.non,data.re))
colnames(cel.file.sample.infor.no.2)=c("filename","subtype")

f.st134 <- factor(cel.file.sample.infor.no.2$subtype)
design.st134 <- model.matrix(~0+f.st134)
colnames(design.st134) <- levels(f.st134)

cancer.data.st134.fpkm <- cbind(rnaSeq.data.non.responding,rnaSeq.data.responding)

cancer.data.st134 <- log2(cbind(rnaSeq.data.non.responding,rnaSeq.data.responding))

fit.st134 <- lmFit(cancer.data.st134, design.st134)
con.trast <- unique(combn(as.character(colnames(design.st134)),2,simplify=FALSE))

con.trast.1 <- unlist(lapply(con.trast, function(u){
  x <- paste(u,collapse = "-")
  x
}))

cont.matrix.st134 <- makeContrasts(contrasts = con.trast.1,levels=design.st134)
fit2.st134  <- contrasts.fit(fit.st134, cont.matrix.st134)
fit2.st134  <- eBayes(fit2.st134,trend=TRUE)

TopTableSt34.all<-topTable(fit2.st134,coef=1,n=dim(fit2.st134)[1])
src.index <- which(rownames(TopTableSt34.all) %in% "SRC")

TopTableSt34.all.1 <- data.frame(Gene=rownames(TopTableSt34.all),TopTableSt34.all)

write.table(TopTableSt34.all.1,file = file.path(output.file.dir,"diff_expr.txt"),
            append = FALSE, quote = F, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = T)


TopTableSt34.all[src.index,]

getDistribution <- function(cancer.data.st134,title,output.file.dir) {
  
  src.index.1 <- which(rownames(cancer.data.st134) %in% "SRC")
  print(cancer.data.st134[src.index.1,])
  X <- as.data.frame(rbind(cbind(log2FPKM=t(cancer.data.st134[src.index.1,1:13]),group=rep("src_non_responding",dim(cancer.data.st134[src.index.1,1:13])[2])),
  cbind(log2FPKM=t(cancer.data.st134[src.index.1,14:28]),group=rep("src_responding",dim(cancer.data.st134[src.index.1,14:28])[2]))))
  X$SRC <- as.numeric(as.character(X$SRC))
  png(file.path(output.file.dir,paste0(title,"_distribution.png")))
  boxplot(X$SRC~X$group,main=title)
  dev.off()
}

getDistribution(cancer.data.st134,title = "log2FPKM_distribution",output.file.dir)
getDistribution(cancer.data.st134.fpkm,title = "FPKM_distribution",output.file.dir)

save.image(file = file.path(output.file.dir,"diff_expr.RData"))
savehistory(file = file.path(output.file.dir,"diff_expr.Rhistory"))