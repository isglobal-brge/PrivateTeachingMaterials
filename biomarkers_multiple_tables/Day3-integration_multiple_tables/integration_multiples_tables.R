## ----setup, echo=FALSE---------------------------------------------------
options(width = 80)
library(knitr)
opts_chunk$set(tidy=FALSE, size='footnotesize', warning=FALSE, cache=TRUE,
               message=FALSE, fig.align='center', out.width='2in')

## ----load_data-----------------------------------------------------------
load("data/breast_TCGA.RData")
group <- droplevels(breast_multi$clin$ER.Status)

## ----data_cc-------------------------------------------------------------
require(CCA)
df1 <- t(breast_multi$RNAseq)[,1:1000]
df2 <- t(breast_multi$RPPA)

## ----cc, eval=FALSE------------------------------------------------------
## resCC <- cc(df1, df2)

## ----rcc, cache=TRUE-----------------------------------------------------
resRCC <- rcc(df1, df2, 0.2, 0.1)

## ----regul_estim, eval=FALSE---------------------------------------------
## regul <- estim.regul(df1, df2)
## resRCC2 <- rcc(df1, df2, regul$lambda1, regul$lambda2)

## ----plot_rcc, fig.show='hide'-------------------------------------------
plt.cc(resRCC)

## ----plot_rcc_out, echo=FALSE--------------------------------------------
plt.cc(resRCC)

## ----multiCCA------------------------------------------------------------
require(PMA)
ddlist <- list(df1, df2)
perm.out <- MultiCCA.permute(ddlist,
                             type=c("standard", "standard"),
                             trace=FALSE)
resMultiCCA <- MultiCCA(ddlist,
                        penalty=perm.out$bestpenalties,
                        ws=perm.out$ws.init,
                        type=c("standard", "standard"),
                        ncomponents=1, trace=FALSE, standardize=TRUE)

## ----multiCCA_out, size='scriptsize'-------------------------------------
rownames(resMultiCCA$ws[[1]]) <- colnames(df1)
rownames(resMultiCCA$ws[[2]]) <- colnames(df2)
head(resMultiCCA$ws[[1]])
head(resMultiCCA$ws[[2]])

## ----getcia--------------------------------------------------------------
library(made4)
library(omicade4)

## ----cia, cache=TRUE-----------------------------------------------------
resCIA <- cia(breast_multi$RNAseq, breast_multi$RPPA)

## ----plot_cia, fig.show='hide'-------------------------------------------
plot(resCIA, classvec=group, nlab=3, clab=0, cpoint=3 )

## ----top_features--------------------------------------------------------
topVar(resCIA, axis=1, topN=5, end="positive")

## ----top_features_neg----------------------------------------------------
topVar(resCIA, axis=1, topN=5, end="negative")

## ----mcia, cache=TRUE----------------------------------------------------
resMCIA <- mcia( breast_multi[ c(1,3,4,5,6,7) ] )

## ----plot_mcia, fig.show='hide'------------------------------------------
plot(resMCIA, axes=1:2, sample.lab=FALSE, sample.legend=FALSE,
     phenovec=group, gene.nlab=2,
     df.color=c("cyan", "magenta", "red4", "brown","yellow", "orange"),
     df.pch=2:7)

## ----top_features_m------------------------------------------------------
topVar(resMCIA, end="positive", axis=1, topN=5)

## ----plot_eigen, fig.show='hide'-----------------------------------------
plot(resMCIA$mcoa$cov2,  xlab = "pseudoeig 1",
     ylab = "pseudoeig 2", pch=19, col="red")
text(resMCIA$mcoa$cov2, labels=rownames(resMCIA$mcoa$cov2),
     cex=1.4, adj=0)

## ----rgcca_pca, eval=FALSE-----------------------------------------------
## # Design matrix C
## # Shrinkage parameters tau = c(tau1, tau2)
## 
## pca.with.rgcca = rgcca(A = list(X, X),
##                        C = matrix(c(0, 1, 1, 0), 2, 2),
##                        tau = c(1, 1))

## ----rgcca_cca, eval=FALSE-----------------------------------------------
## # X1 = Block1 and X2 = Block2
## # Design matrix C
## # Shrinkage parameters tau = c(tau1, tau2)
## 
## cca.with.rgcca = rgcca(A= list(X1, X2),
##                        C = matrix(c(0, 1, 1, 0), 2, 2),
##                        tau = c(0, 0))

## ----rgcca_pls, eval=FALSE-----------------------------------------------
## # X1 = Block1 and X2 = Block2
## # Design matrix C
## # Shrinkage parameters tau = c(tau1, tau2)
## 
## pls.with.rgcca = rgcca(A= list(X1, X2),
##                        C = matrix(c(0, 1, 1, 0), 2, 2),
##                        tau = c(1, 1))

## ----rgcca_rda, eval=FALSE-----------------------------------------------
## # X1 = Block1 and X2 = Block2
## # Design matrix C
## # Shrinkage parameters tau = c(tau1, tau2)
## 
## ra.with.rgcca = rgcca(A= list(X1, X2),
##                        C = matrix(c(0, 1, 1, 0), 2, 2),
##                        tau = c(1, 0))

## ----gcca_rgcca, eval = FALSE--------------------------------------------
## # X1 = Block1, ..., XJ = BlockJ, X_{J+1} = [X1, ..., XJ]
## # (J+1)*(J+1) Design matrix C
## C = matrix(c(0, 0, 0, ..., 0, 1,
##              0, 0, 0, ..., 0, 1,
##                     ...
##              1, 1, 1, ..., 1, 0), J+1, J+1)
## # Shrinkage parameters tau = c(tau1, ...,  tauJ, tau_{J+1})
## gcca.with.rgcca = rgcca(A= list(X1, ..., XJ, cbind(X1, ..., XJ)),
##                        C = C, tau = rep(0, J+1),
##                        scheme = "factorial")

## ----mcia_rgcca, eval = FALSE--------------------------------------------
## # X1 = Block1, ..., XJ = BlockJ, X_{J+1} = [X1, ..., XJ]
## # (J+1)*(J+1) Design matrix C
## C = matrix(c(0, 0, 0, ..., 0, 1,
##              0, 0, 0, ..., 0, 1,
##                   ...
##              1, 1, 1, ..., 1, 0), J+1, J+1)
## # Shrinkage parameters tau = c(tau1, ...,  tauJ, tau_{J+1})
## mcoa.with.rgcca = rgcca(A= list(X1, ..., XJ, cbind(X1, ..., XJ)),
##                        C = C, tau = c(rep(1, J), 0),
##                        scheme = "factorial")

## ----rgcca_1-------------------------------------------------------------
library(RGCCA)
load("data/breast_TCGA.RData")
X <- t(breast_multi$RNAseq)
Y <- t(breast_multi$miRNA)
Z <- t(breast_multi$RPPA)
A <- list(rnaseq=X, miRNA=Y, RPPA=Z)
A <- lapply(A, scale)

## ----C-------------------------------------------------------------------
C <- matrix(c(1,0,1,0,1,1,0,0,1), nrow=3, byrow = TRUE)
C

## ----rgcca.factorial, cache=TRUE-----------------------------------------
rgcca.factorial <- rgcca(A, C=C, tau = rep(0, 3), 
                           scheme ="factorial", ncomp=c(2,2,2),
                           scale = FALSE, verbose = FALSE)

## ----sol, eval=FALSE-----------------------------------------------------
## rgcca.factorial$a # weight vectors

## ----block---------------------------------------------------------------
Y.block <- rgcca.factorial$Y
lapply(Y.block, head)

## ----ave-----------------------------------------------------------------
rgcca.factorial$AVE

## ----rgcca_ind, fig.show='hide'------------------------------------------
source("R/plotInd.R")
plotInd(rgcca.factorial, group)

## ----miRNA---------------------------------------------------------------
source("R/topVars.R")
topVars(rgcca.factorial, axis=1, end="pos", topN=5)

## ----prev, cache=TRUE----------------------------------------------------
X <- t(breast_multi$RNAseq)
Y <- t(breast_multi$miRNA)
AA <- list(rnaseq=X, miRNA=Y, group=as.numeric(group)-1) 
C <-  matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
C

## ----sgcca, cache=TRUE---------------------------------------------------
sgcca.breast = sgcca(AA, C, c1 = c(.071,.2, 1),
                     ncomp = c(1, 1, 1),
                     scheme = "centroid",
                     scale = TRUE,
                     verbose = FALSE)

## ----sgcca1--------------------------------------------------------------
ss1 <- sgcca.breast$a[[1]] != 0
sum(ss1)
features.miRNA <- rownames(sgcca.breast$a[[1]])[ss1]
head(features.miRNA)

## ----sgcca2--------------------------------------------------------------
ss2 <- sgcca.breast$a[[2]] != 0
sum(ss2)
features.miRNA <- rownames(sgcca.breast$a[[2]])[ss2]
head(features.miRNA)

## ----sgcca3--------------------------------------------------------------
source("R/selectVars.R")
ans <- selectVars(sgcca.breast, table=1, end = "pos")
length(ans)
ans

## ----sgcca4--------------------------------------------------------------
ans <- selectVars(sgcca.breast, table=2, end = "neg")
length(ans)
ans

