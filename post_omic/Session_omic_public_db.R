#
# GEO: predict cronological age using DNA methylation
#

library(GEOquery)
library(glmnet)

gse58045 <- getGEO("GSE58045")[[1]]

setwd("c:/juan/CREAL/GitHub/PrivateTeachingMaterials/post_omic/")
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


# external validation
gse41037 <- getGEO("GSE41037")[[1]]

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
library(TCGAutils)
library(omicade4)


###  Breast

# download data
brca <- curatedTCGAData("BRCA", c("miRNASeqGene"), FALSE)
sampleTables(brca)

# select tumor samples
brca.c <- splitAssays(brca, c("01"))

# prepare survival data
pheno <- colData(brca.c)[, getClinicalNames("BRCA")]
pheno$time <- apply(pheno[, c("days_to_death", "days_to_last_followup")],
                    1, function(x) max(x, na.rm=TRUE))
pheno$time[pheno$time<0] <- 0
pheno$status <- pheno$vital_status

# get miRNA
mirna <- assay(brca.c)
colnames(mirna) <- substr(colnames(mirna),1,12)


# create data.frame with complete cases
o <- intersect(colnames(mirna), rownames(pheno))
pheno.o <- pheno[o, ]
mirna.o <- mirna[, o]
breast.df <- data.frame(pheno.o, t(mirna.o))
load("data/breast.Rdata")


# run survival analysis
library(survival)
ii <- grep("hsa", names(breast.df))

ff <- function(i, data) {
  mod <- coxph(Surv(time, status) ~ data[,i], data=data)
  ans <- summary(mod)$coefficients[5]
  ans
}

out <- sapply(ii, ff, data=breast.df)

res <- data.frame(miRNA=names(breast.df)[ii],
                  pval=out)
res$padj <- p.adjust(res$pval, method="BH")
head(res[order(res$pval),])



# gliobastoma


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

