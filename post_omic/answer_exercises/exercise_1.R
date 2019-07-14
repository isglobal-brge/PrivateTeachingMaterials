#
# EXERCISE 1
#

library(Biobase)
library(SummarizedExperiment)
load(file="data_exercises/breast_tcga.RData")
breast

data <- colData(breast)
table(breast$er)

boxplot(assays(breast)$counts[1,] ~ breast$er)


roi <- GRanges(seqnames="chr6", IRanges(start=151.2e6,
                                 end=151.8e6))
breast.roi <- subsetByOverlaps(breast, roi)
assay(breast.roi)
rowRanges(breast.roi)
