#
# GEO: predict cronological age using DNA methylation
#

library(GEOquery)
library(glmnet)


setwd("c:/Juan/CREAL/BayesianPrediction/Bayesian_clock/paper")
load("data/GSE58045.Rdata")
gse58045
names(pData(gse58045))

# Outcome variable
y <- as.numeric(pData(gse58045)$`age:ch1`)

# Predictor variables
x <- t(exprs(gse58045))
x.imp <- impute::impute.knn(x)$data

# elastic net

cv <- cv.glmnet(x.imp, y, alpha = 0.5)
cv$lambda.min

model <- glmnet(x.imp, y, alpha = 0.5, 
                      lambda = cv$lambda.min)
coefs <- coef(model)[,1]
coefMod <- coefs[coefs!=0]
CpGs <- names(coefs)[coefs!=0]

load("data/GSE41037.Rdata")
common <- intersect(CpGs, rownames(gse41037))
validation <- t(exprs(gse41037))[, common]
age.pred <- coefMod[1] + validation%*%coefMod[common]

age.validation <- as.numeric(pData(gse41037)$`age:ch1`)
mod <- lm(age.validation ~ age.pred[,1])
plot(age.pred[,1], age.validation, xlab="DNAmAge", ylab="Chronological Age")
abline(mod, col="blue")


#
# TCGA
#


library(curatedTCGAData)
library(MultiAssayExperiment)
library(omicade4)

gbm <- curatedTCGAData("GBM", c("RNASeq2GeneNorm", 
                                   "RPPAArray", "miRNAArray"), FALSE)

gbmL <- assays(gbm)
lapply(gbmL, dim)
gbmL <- lapply(gbmL, function(m){
  m <- data.matrix(m)
  colnames(m) <- substring(colnames(m), 1, 12)
  ## Remove features with missings
  m[ rowSums(is.na(m)) == 0, ]
})

## Log transform count data
gbmL[[2]] <- log2(gbmL[[2]] + 1)

## select common cases
common <- intersect(colnames(gbmL[[1]]), colnames(gbmL[[2]]))
common <- intersect(common, colnames(gbmL[[3]]))
length(common)

# remove non-expressed genes
mask <- !apply(gbmL[[2]][,common], 1, function(x) all(x==0))

gbmL.common <- list(miRNA = gbmL[[1]][,common],
                    RNAseq = gbmL[[2]][mask,common],
                    RPPA = gbmL[[3]][,common])

load("data/gbm.Rdata")

res <- mcia(gbmL.common)
plot(res, axes=1:2, sample.lab=FALSE, sample.legend=FALSE, gene.nlab=2,
     df.color=c("cyan", "magenta", "red4"),
     df.pch=2:4)

