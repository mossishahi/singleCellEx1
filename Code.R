library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)

data <- read.delim(file = "PBMC_scRNA-seq.txt", sep = " ")
object <- SingleCellExperiment(assays = list(counts = as.matrix(data)))
isSpike(object, "ERCC") <- grepl("ERC", rownames(object))
isSpike(object, "MT") <- grepl("MT-", rownames(object))
QC <- calculateQCMetrics(object , feature_controls = list(ERCC = isSpike(object, "ERCC"), MT = isSpike(object, "MT")))
cells <- colData(QC)
cellsdata <- data.frame(cells@listData)

ggplot(cellsdata, aes(x = cellsdata$total_counts)) + geom_histogram(binwidth = 20)
ggplot(cellsdata, aes(x = cellsdata$total_features)) + geom_histogram(binwidth = 20)
ggplot(data.frame(cellsdata$pct_counts_ERCC[cellsdata$pct_counts_ERCC != 0]), aes(x = cellsdata$pct_counts_ERCC[cellsdata$pct_counts_ERCC != 0])) + geom_histogram(binwidth = 0.001)
ggplot(cellsdata, aes(x = cellsdata$pct_counts_MT)) + geom_histogram(binwidth = 0.01)

filter1 <- (cellsdata$total_counts > 3*mad(cellsdata$total_counts))
filter2 <- (cellsdata$total_features > 3*mad(cellsdata$total_features))
filter3 <- (cellsdata$pct_counts_ERCC < 3*mad(cellsdata$pct_counts_ERCC[cellsdata$pct_counts_ERCC != 0]))
filter4 <- (cellsdata$pct_counts_MT < 3*mad(cellsdata$pct_counts_MT))

newdatacells <- data[,filter1 & filter2 & filter3 & filter4]

avgexp <- calcAverage(object)
ggplot(data.frame(avgexp), aes(x = avgexp)) + geom_histogram(binwidth = 0.01)

genefilter <- (avgexp > 0.07)
newdatagenes <- newdatacells[genefilter,]

newobject <- SingleCellExperiment(assays = list(counts = as.matrix(newdatagenes)))
isSpike(newobject, "ERCC") <- grepl("ERC", rownames(newobject))
isSpike(newobject, "MT") <- grepl("MT-", rownames(newobject))
newQC <- calculateQCMetrics(newobject , feature_controls = list(ERCC = isSpike(newobject, "ERCC"), MT = isSpike(newobject, "MT")))

qclust <- quickCluster(newQC, min.size = 30)
normalQC <- computeSumFactors(newQC, sizes = 15, clusters = qclust)
normalQC <- normalize(normalQC)

spikenormalQC <- computeSpikeFactors(newQC, type = c("MT", "ERCC"))
spikenormalQC <- normalize(spikenormalQC)

plotPCA(
  spikenormalQC,
  exprs_values = "counts",
  size_by = "total_features"
)

library(Rtsne)
plotTSNE(
  spikenormalQC,
  exprs_values = "counts",
  perplexity = 130,
  size_by = "total_features",
  rand_seed = 123456
)

var.fit <- trendVar(spikenormalQC, parametric=TRUE, use.spikes=FALSE, span=0.2)
var.out <- decomposeVar(spikenormalQC, var.fit)
cur.spike <- isSpike(spikenormalQC)
spikenormalQC <- denoisePCA(spikenormalQC, technical=var.fit$trend, min.rank = 2)
pcs <- reducedDim(spikenormalQC)

idx <- kmeans(pcs,centers = 5,iter.max = 100 , nstart = 200, algorithm = c("Hartigan-Wong", "Lloyd", "Forgy","MacQueen"), trace=FALSE)
idx = factor(idx$cluster)
spikenormalQC$cluster5 <- idx
plotReducedDim(spikenormalQC, use_dimred="PCA" ,colour_by="cluster5" )
