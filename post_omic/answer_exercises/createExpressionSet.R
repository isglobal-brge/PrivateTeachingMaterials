
library(Biobase)
id <- paste0("id", 1:10)
data <- data.frame(id=id, age=1:10, row.names=id)
genes <- matrix(rnorm(100), nrow=10)
colnames(genes) <- id
phenoData <- AnnotatedDataFrame(data=data)
exprs <- ExpressionSet(genes, phenoData=phenoData,
                       annotation="")
exprs
