### R code from vignette source '/Library/Frameworks/R.framework/Versions/3.5/Resources/library/HTqPCR/doc/HTqPCR.Rnw'

###################################################
### code chunk number 1: Prepare parameters
###################################################
options(width=65)
set.seed(123)


###################################################
### code chunk number 2: Load package
###################################################
library("HTqPCR")


###################################################
### code chunk number 3: Extract R code (eval = FALSE)
###################################################
## all.R.commands <- system.file("doc", "HTqPCR.Rnw", package = "HTqPCR")
## Stangle(all.R.commands)


###################################################
### code chunk number 4: All functions
###################################################
ls("package:HTqPCR")


###################################################
### code chunk number 5: Load example data
###################################################
data(qPCRraw)
data(qPCRpros)
class(qPCRraw)


###################################################
### code chunk number 6: Information contained in qPCRsets
###################################################
slotNames(qPCRraw)
phenoData(qPCRraw)
pData(qPCRraw)
pData(qPCRraw) <- data.frame(Genotype=rep(c("A", "B"), each=3), Replicate=rep(1:3, 2))
pData(qPCRraw)
featureData(qPCRraw)
head(fData(qPCRraw))


###################################################
### code chunk number 7: Example input files
###################################################
path <- system.file("exData", package="HTqPCR")
head(read.delim(file.path(path, "files.txt")))


###################################################
### code chunk number 8: Read raw data
###################################################
files <- read.delim(file.path(path, "files.txt"))
raw <- readCtData(files=files$File, path=path)


###################################################
### code chunk number 9: Show qPCRset data object
###################################################
show(raw)


###################################################
### code chunk number 10: Ct overview ex 1
###################################################
g <- featureNames(raw)[1:10]
plotCtOverview(raw, genes=g, xlim=c(0,50), groups=files$Treatment, conf.int=TRUE,
ylim=c(0,55))


###################################################
### code chunk number 11: Ct overview ex 2
###################################################
plotCtOverview(raw, genes=g, xlim=c(0,50), groups=files$Treatment,
calibrator="Control")


###################################################
### code chunk number 12: HTqPCR.Rnw:242-245
###################################################
par(mfrow=c(2,1))
g <- featureNames(raw)[1:10]
plotCtOverview(raw, genes=g, xlim=c(0,50), groups=files$Treatment, conf.int=TRUE,
ylim=c(0,55))
plotCtOverview(raw, genes=g, xlim=c(0,50), groups=files$Treatment,
calibrator="Control")


###################################################
### code chunk number 13: Ct card ex 1
###################################################
plotCtCard(raw, col.range=c(10,35), well.size=2.6)


###################################################
### code chunk number 14: Ct card ex 2
###################################################
featureClass(raw) <- factor(c("Marker", "TF", "Kinase")[sample(c(1,1,2,2,1,3),
384, replace=TRUE)])
plotCtCard(raw, plot="class", well.size=2.6)


###################################################
### code chunk number 15: HTqPCR.Rnw:267-268
###################################################
plotCtCard(raw, col.range=c(10,35), well.size=2.6)


###################################################
### code chunk number 16: HTqPCR.Rnw:270-271
###################################################
featureClass(raw) <- factor(c("Marker", "TF", "Kinase")[sample(c(1,1,2,2,1,3),
384, replace=TRUE)])
plotCtCard(raw, plot="class", well.size=2.6)


###################################################
### code chunk number 17: Ct replicates
###################################################
plotCtReps(qPCRraw, card=2, percent=20)


###################################################
### code chunk number 18: HTqPCR.Rnw:288-289
###################################################
plotCtReps(qPCRraw, card=2, percent=20)


###################################################
### code chunk number 19: Ct variation ex 1
###################################################
raw.mix <- raw
#exprs(raw.mix)[,6] <- sample(exprs(raw[,6]))
plotCtVariation(raw.mix, variation="sd", log=TRUE, main="SD of replicated features", col="lightgrey")


###################################################
### code chunk number 20: Ct variation ex 2
###################################################
raw.variation <- plotCtVariation(raw.mix, type="detail", add.featurenames=TRUE, pch=" ", cex=1.2)


###################################################
### code chunk number 21: Ct variation ex 2 in detail
###################################################
names(raw.variation)
head(raw.variation[["Var"]][,1:4])
head(raw.variation[["Mean"]][,1:4])
apply(raw.variation[["Var"]][,3:7], 2, summary)
colSums(raw.variation[["Var"]][,3:7]>20)


###################################################
### code chunk number 22: HTqPCR.Rnw:328-329
###################################################
raw.mix <- raw
#exprs(raw.mix)[,6] <- sample(exprs(raw[,6]))
plotCtVariation(raw.mix, variation="sd", log=TRUE, main="SD of replicated features", col="lightgrey")


###################################################
### code chunk number 23: HTqPCR.Rnw:331-332
###################################################
raw.variation <- plotCtVariation(raw.mix, type="detail", add.featurenames=TRUE, pch=" ", cex=1.2)


###################################################
### code chunk number 24: Plot Ct categories ex 1
###################################################
raw.cat <- raw
plotCtCategory(raw.cat)


###################################################
### code chunk number 25: Plot Ct categories ex 2
###################################################
plotCtCategory(raw.cat, stratify="class")


###################################################
### code chunk number 26: HTqPCR.Rnw:373-376
###################################################
par(mfrow=c(2,1))
raw.cat <- raw
plotCtCategory(raw.cat)
plotCtCategory(raw.cat, stratify="class")


###################################################
### code chunk number 27: Plot Ct categories ex 3
###################################################
plotCtCategory(raw.cat, by.feature=TRUE, cexRow=0.1)


###################################################
### code chunk number 28: HTqPCR.Rnw:391-392
###################################################
plotCtCategory(raw.cat, by.feature=TRUE, cexRow=0.1)


###################################################
### code chunk number 29: Normalise data
###################################################
q.norm <- normalizeCtData(raw.cat, norm="quantile")
sr.norm <- normalizeCtData(raw.cat, norm="scale.rank")
nr.norm <- normalizeCtData(raw.cat, norm="norm.rank")
d.norm <- normalizeCtData(raw.cat, norm="deltaCt", deltaCt.genes=c("Gene1", "Gene60"))
g.norm <- normalizeCtData(raw.cat, norm="geometric.mean")


###################################################
### code chunk number 30: Normalisation comparison
###################################################
plot(exprs(raw), exprs(q.norm), pch=20, main="Quantile normalisation", col=rep(brewer.pal(6, "Spectral"), each=384))


###################################################
### code chunk number 31: HTqPCR.Rnw:444-460
###################################################
col <- rep(brewer.pal(6, "Spectral"), each=384)
col2 <- brewer.pal(5, "Dark2")
par(mfrow=c(3,2), mar=c(2,2,2,2))
# All methods individually
plot(exprs(raw), exprs(q.norm), pch=20, main="Quantile normalisation", col=col)
plot(exprs(raw), exprs(sr.norm), pch=20, main="Rank invariant scaling", col=col)
plot(exprs(raw), exprs(nr.norm), pch=20, main="Rank invariant normalisation", col=col)
plot(exprs(raw), exprs(d.norm), pch=20, main="deltaCt normalisation", col=col)
plot(exprs(raw), exprs(g.norm), pch=20, main="Geometric mean normalisation", col=col)
# Just a single sample, across methods
plot(exprs(raw)[,3], exprs(q.norm)[,3], pch=20, col=col2[1], main="Comparison of methods for sample 3", ylim=c(-10,40))
points(exprs(raw)[,3], exprs(sr.norm)[,3], pch=20, col=col2[2])
points(exprs(raw)[,3], exprs(nr.norm)[,3], pch=20, col=col2[3])
points(exprs(raw)[,3], exprs(d.norm)[,3], pch=20, col=col2[4])
points(exprs(raw)[,3], exprs(g.norm)[,3], pch=20, col=col2[5])
legend(8, 40, legend=c("Quantile", "Rank.invariant scaling", "Rank.invariant normalization", "deltaCt", "Geometric.mean"), col=col2, lwd=2, bty="n")


###################################################
### code chunk number 32: Ct correlations
###################################################
plotCtCor(raw, main="Ct correlation")


###################################################
### code chunk number 33: HTqPCR.Rnw:540-541
###################################################
plotCtCor(raw, main="Ct correlation")


###################################################
### code chunk number 34: Summary of Ct values
###################################################
summary(raw)


###################################################
### code chunk number 35: Ct density
###################################################
plotCtDensity(sr.norm)


###################################################
### code chunk number 36: Ct histogram
###################################################
plotCtHistogram(sr.norm)


###################################################
### code chunk number 37: HTqPCR.Rnw:571-574
###################################################
par(mfrow=c(1,2), mar=c(3,3,2,1))
plotCtDensity(sr.norm)
plotCtHistogram(sr.norm)


###################################################
### code chunk number 38: HTqPCR.Rnw:585-592
###################################################
par(mfrow=c(3,2), mar=c(2,2,2,1))
plotCtDensity(qPCRraw, main="Raw Ct values")
plotCtDensity(q.norm, main="quantile")
plotCtDensity(sr.norm, main="scale.rankinvariant")
plotCtDensity(nr.norm, main="norm.rankinvariant")
plotCtDensity(d.norm, main="deltaCt")
plotCtDensity(g.norm, main="geometric.mean")


###################################################
### code chunk number 39: Ct boxes
###################################################
plotCtBoxes(sr.norm, stratify="class")


###################################################
### code chunk number 40: HTqPCR.Rnw:607-608
###################################################
plotCtBoxes(sr.norm, stratify="class")


###################################################
### code chunk number 41: Ct scatter ex 1
###################################################
plotCtScatter(sr.norm, cards=c(1,2), col="type", diag=TRUE)


###################################################
### code chunk number 42: Ct scatter ex 2
###################################################
plotCtScatter(sr.norm, cards=c(1,4), col="class", diag=TRUE)


###################################################
### code chunk number 43: HTqPCR.Rnw:628-631
###################################################
par(mfrow=c(1,2), mar=c(3,3,2,1))
plotCtScatter(sr.norm, cards=c(1,2), col="type", diag=TRUE)
plotCtScatter(sr.norm, cards=c(1,4), col="class", diag=TRUE)


###################################################
### code chunk number 44: Ct pairs
###################################################
plotCtPairs(sr.norm, col="type", diag=TRUE)


###################################################
### code chunk number 45: HTqPCR.Rnw:649-650
###################################################
plotCtPairs(sr.norm, col="type", diag=TRUE)


###################################################
### code chunk number 46: Ct heatmap
###################################################
plotCtHeatmap(raw, gene.names="", dist="euclidean")


###################################################
### code chunk number 47: HTqPCR.Rnw:669-670
###################################################
plotCtHeatmap(raw, gene.names="", dist="euclidean")


###################################################
### code chunk number 48: CV across samples
###################################################
plotCVBoxes(qPCRraw, stratify="class")
plotCVBoxes(qPCRraw, stratify="type")


###################################################
### code chunk number 49: HTqPCR.Rnw:688-691
###################################################
par(mfrow=c(1,2), mar=c(2,2,2,1))
plotCVBoxes(qPCRraw, stratify="class")
plotCVBoxes(qPCRraw, stratify="type")


###################################################
### code chunk number 50: Cluster Ct
###################################################
clusterCt(sr.norm, type="samples")


###################################################
### code chunk number 51: HTqPCR.Rnw:717-718
###################################################
clusterCt(sr.norm, type="samples")


###################################################
### code chunk number 52: Plot subclusters
###################################################
cluster.list <- clusterCt(sr.norm, type="genes", n.cluster=6, cex=0.5)


###################################################
### code chunk number 53: HTqPCR.Rnw:737-738
###################################################
cluster.list <- clusterCt(sr.norm, type="genes", n.cluster=6, cex=0.5)


###################################################
### code chunk number 54: Principal components analysis
###################################################
plotCtPCA(qPCRraw)
plotCtPCA(qPCRraw, features=FALSE)


###################################################
### code chunk number 55: HTqPCR.Rnw:756-759
###################################################
par(mfrow=c(1,2), mar=c(2,2,2,1))
plotCtPCA(qPCRraw)
plotCtPCA(qPCRraw, features=FALSE)


###################################################
### code chunk number 56: Object history
###################################################
getCtHistory(sr.norm)


###################################################
### code chunk number 57: Multiple samples per card
###################################################
# Example with 2 or 4 samples per 384 well card.
sample2.order <- rep(c("subSampleA", "subSampleB"), each=192)
sample4.order <- rep(c("subA", "subB", "subC", "subD"), each=96)
# Splitting the data into all individual samples
qPCRnew2 <- changeCtLayout(sr.norm, sample.order=sample2.order)
show(qPCRnew2)
qPCRnew4 <- changeCtLayout(sr.norm, sample.order=sample4.order)
show(qPCRnew4)


###################################################
### code chunk number 58: Card history
###################################################
getCtHistory(qPCRnew4)


###################################################
### code chunk number 59: Example SDS data
###################################################
path <- system.file("exData", package="HTqPCR")
cat(paste(readLines(file.path(path, "SDS_sample.txt"), n=19), "\n"))


###################################################
### code chunk number 60: Example SDS data 2
###################################################
readLines(file.path(path, "SDS_sample.txt"), n=20)


###################################################
### code chunk number 61: SDS format
###################################################
path <- system.file("exData", package = "HTqPCR")
raw <- readCtData(files="SDS_sample.txt", path=path, format="SDS")
show(raw)


###################################################
### code chunk number 62: SDS format
###################################################
path <- system.file("exData", package = "HTqPCR")
raw <- readCtData(files="LightCycler_sample.txt", path=path, format="LightCycler")
show(raw)


###################################################
### code chunk number 63: CFX format
###################################################
path <- system.file("exData", package = "HTqPCR")
raw <- readCtData(files="CFX_sample.txt", path=path, format="CFX", n.features=330)
show(raw)


###################################################
### code chunk number 64: BioMark format
###################################################
exPath <- system.file("exData", package="HTqPCR")
raw1 <- readCtData(files="BioMark_sample.csv", path=exPath, format="BioMark", n.features=48, n.data=48)
dim(raw1)
raw2 <- readCtData(files="BioMark_sample.csv", path=exPath, format="BioMark", n.features=48*48, n.data=1)
dim(raw2)


###################################################
### code chunk number 65: Ct microfluidic array
###################################################
plotCtArray(raw1)


###################################################
### code chunk number 66: HTqPCR.Rnw:1078-1079
###################################################
plotCtArray(raw1)


###################################################
### code chunk number 67: OpenArray format
###################################################
exPath <- system.file("exData", package="HTqPCR")
raw1 <- readCtData(files="OpenArray_sample.csv", path=exPath, format="OpenArray", n.features=846, n.data=6)
dim(raw1)
raw2 <- readCtData(files="OpenArray_sample.csv", path=exPath, format="OpenArray", n.features=846*6, n.data=1)
dim(raw2)


###################################################
### code chunk number 68: Check HTqPCR news
###################################################
news(Version>1.7, package="HTqPCR")


###################################################
### code chunk number 69: sessionInfo
###################################################
toLatex(sessionInfo())


