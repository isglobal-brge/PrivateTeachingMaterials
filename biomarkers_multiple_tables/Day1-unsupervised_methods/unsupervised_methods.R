## ----data----------------------------------------------------------------
require(graphics)
data(USArrests)
head(USArrests)

## ----pca-----------------------------------------------------------------
princomp(USArrests)
mod <- princomp(USArrests)
summary(mod)

## ----pca2----------------------------------------------------------------
mod2 <- princomp(USArrests)
summary(mod2)

## ----plotPCA-------------------------------------------------------------
plot(mod2)

## ----pca3----------------------------------------------------------------
loadings(mod2)

## ----pcaPlot3------------------------------------------------------------
biplot(mod2)

## ----plotPCA2------------------------------------------------------------
library(made4)
load("data/breast_TCGA.RData")
group <- droplevels(breast_multi$clin$ER.Status)
rnaseq <- breast_multi$RNAseq
out <- ord(rnaseq, trans=FALSE, type="pca", classvec=group)
plot(out, nlab=3, arraylabels=rep("T", 79))

## ----plotPCA2b-----------------------------------------------------------
par(mfrow=c(2,1))
plotarrays(out$ord$co, classvec=group)
plotgenes(out, col="blue")

## ----list_genes----------------------------------------------------------
ax1 <- topgenes(out, axis=1, n=5, ends="pos")
ax2 <- topgenes(out, axis=2, n=5, ends="neg")
cbind(pos=ax1, neg=ax2)

## ----scale---------------------------------------------------------------
  rnaseq.s <- scale(rnaseq, center= TRUE, scale = TRUE)

## ----impute--------------------------------------------------------------
library(impute)
rnaseq.s.imp <- impute.knn(rnaseq.s, rowmax = 0.5,
                         colmax = 0.8)$data #samples in columns!!!

## ----plotPCA3------------------------------------------------------------
out <- ord(rnaseq.s.imp, trans=FALSE, type="pca", classvec=group)
plot(out, nlab=3, arraylabels=rep("T", 79))

## ----var_exp-------------------------------------------------------------
summary(out$ord)

## ----boot, fig.show='hide'-----------------------------------------------
library(nFactors)
ev <- eigen(cor(USArrests)) # get eigenvalues
ap <- parallel(subject=nrow(USArrests), var=ncol(USArrests),
  rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS) 

## ----lcpm----------------------------------------------------------------
require(LPCM)
data(gaia)
dim(gaia)
names(gaia)

## ----lcpmPlot------------------------------------------------------------
 lpc1 <-  lpc(gaia[,7:8])
 plot(lpc1, curvecol="red", lwd=3)

## ----lcpmPlot2-----------------------------------------------------------
require(scatterplot3d)
lpc2 <- lpc(gaia[,7:9])
plot(lpc2, curvecol=2, type=c("curve","mass"))

## ------------------------------------------------------------------------
dd <- read.delim("data/cars.tab")
head(dd)

## ----norm----------------------------------------------------------------
dd.ok <- dd[ , -c(1,2)]
medians <- apply(dd.ok, 2, median)
mads <- apply(dd.ok, 2, mad) # median absolute deviation
dd.ok <- scale(dd.ok, center=medians, scale=mads) 

## ----dist----------------------------------------------------------------
dd.dist <- dist(dd.ok)
dd.dist.camb <- dist(dd.ok, method="canberra")

## ----dist2---------------------------------------------------------------
library(cluster)
dd.dist.gower <- daisy(dd.ok, metric="gower") 

## ----hclust--------------------------------------------------------------
dd.hclust <- hclust(dd.dist)
plot(dd.hclust, labels=dd$Car, main="")

## ----hclustCut-----------------------------------------------------------
groups3.hclust <- as.factor(cutree(dd.hclust, 3))
table(groups3.hclust)

## ----hclustGroup---------------------------------------------------------
dd$Car[groups3.hclust==1]

## ----hclustGroup2--------------------------------------------------------
sapply(unique(groups3.hclust), function (x) dd$Car[groups3.hclust==x])

## ----hclustCol-----------------------------------------------------------
dend <- as.dendrogram(dd.hclust)
dend2 <- dendextend::color_labels(dend, k=3)
dendextend::labels(dend2) <- dd$Car
plot(dend2)

## ----pam-----------------------------------------------------------------
require(cluster)
dd.pam <- pam(dd.dist, 3)
dd.kmeans <- kmeans(dd.dist, 3)
groups3.pam <- as.factor(dd.pam$clustering)
groups3.kmeans <-  as.factor(dd.kmeans$cluster)

## ----comparison----------------------------------------------------------
table(groups3.pam, groups3.hclust)
table(groups3.pam, groups3.kmeans)

## ----pamPlot-------------------------------------------------------------
plot(dd.pam)

## ----pamPlotAny----------------------------------------------------------
plot(silhouette(cutree(dd.hclust,4), dd.dist))

## ----mclust--------------------------------------------------------------
require(mclust)
dd.mclust <- Mclust(dd.ok)
summary(dd.mclust)

## ----mclustG-------------------------------------------------------------
groups3.mclust <- Mclust(dd.ok, G=3)$class
table(groups3.mclust, groups3.pam)

## ----mclustFig2D---------------------------------------------------------
mclust2Dplot(dd.ok[,1:2], parameters=dd.mclust$parameters, 
             z=dd.mclust$z, what = "classification")

## ----selectK-------------------------------------------------------------
require(vegan)
k.cal <- cascadeKM(dd.ok, inf.gr=2, sup.gr=6, criterion="calinski")
k.ssi <- cascadeKM(dd.ok, inf.gr=2, sup.gr=6, criterion="ssi")

## ----plotCal-------------------------------------------------------------
plot(k.cal)

## ----plotSsi-------------------------------------------------------------
plot(k.ssi)

## ----plotTest------------------------------------------------------------
require(sigclust)
mod.sig <-sigclust(dd.ok, nsim=1000)
plot(mod.sig, arg="pvalue")

## ----validation----------------------------------------------------------
require(clValid)
val.intern <- clValid(dd.ok, 2:6, clMethods = c("hierarchical",
 "kmeans", "diana", "pam",  "model"), validation = "internal")
optimalScores(val.intern)

## ----big-----------------------------------------------------------------
require(fastcluster)
hclust.fast <- hclust(dd.dist)

