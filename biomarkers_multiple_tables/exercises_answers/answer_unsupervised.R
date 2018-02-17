library(made4)

#
# task 1
#

setwd("c:/Juan/CREAL/GitHub/biomarkers_multiple_tables")
load("data_exercises/nhanes.Rdata")
ls()

# select variables
dd <- nhanes.nut[,1:29]
group <- as.factor(nhanes.nut$CHOL2)
table(group)


# scale data
dd.s <- lapply(dd, scale)

#coa
mod <- ord(dd, type="coa", trans=TRUE, classvec = group)
plot(mod, nlab=5)

# top variables second axis (because it seems that fist axis is not separating
# normal vs hypercol
topgenes(mod, axis=2, end="pos")
topgenes(mod, axis=2, end="neg")

# Let us investigate whether there is association between second axis and CHOL2 (variable group)
pca1 <- out$ord$co[,2]
mod <- glm(group ~ pca1, family="binomial")
summary(mod)


#
# task 2
#

library(made4)
setwd("c:/Juan/CREAL/GitHub/biomarkers_multiple_tables")
load("data_exercises/nci60.Rdata")
miRNA <- nci60$miRNA

pca.miRNA <- ord(miRNA, classvec=cancer)
plot(pca.miRNA, nlab=5)

# these are associated with Leukemia 
# (notice in the figure that in this case leukemias are in the negative side)
topgenes(pca.miRNA, axis=1, n=5, ends="neg")

# these are associated with Melanoma
# (notice in the figure that in this case melanomas are in the negative side)
ax1.pos <- topgenes(pca.miRNA, axis=1, n=5, ends="pos")
ax1.pos
ax2.pos <- topgenes(pca.miRNA, axis=2, n=5, ends="pos")
intersect(ax1.pos, ax2.pos)

# these are associated with Central Nervous System
topgenes(pca.miRNA, axis=2, n=5, ends="neg")

# variability
summary(pca.miRNA$ord)


plot(pca.miRNA, axis1 = 3, axis2=4)

#
# task 3
#

library(cluster)
setwd("c:/Juan/CREAL/GitHub/biomarkers_multiple_tables")
load("data_exercises/diet.Rdata")


# This is another question raised after some discussion
# Perform clustering on variables!!! 
X <- cbind(X2, X3)
# X <- X[complete.cases(X), ]
X.s <- scale(X)
dd <- dist(t(X.s))
dd.h <- hclust(dd)
plot(dd.h)

# Compare what happens if data are not standardized
X <- cbind(X2, X3)
dd2 <- dist(t(X))
dd2.h <- hclust(dd2)
plot(dd2.h)
summary(X$t_energy)
summary(X$t_folic)

# Assign individuals to each group and perform association study
library(Hmisc)
gg <- cutree(dd.h, 2)
gg
dim(X)
cut2(X$t_folic,3)
ff <- function(x){
  as.numeric(as.factor(cut2(x, g=4)))
}
X.q <- apply(X, 2, ff)
cbind(X[1:5,1:3], X.q[1:5,1:3])

X.q.nomed <- X.q[, gg==1]
X.q.med <- X.q[, gg==2]
score.nomed <- apply(X.q.nomed, 1, sum)
score.med <- apply(X.q.med, 1, sum)

mod <- glm(casoc ~ score.med, data=X1, family="binomial")
summary(mod)
1-exp(coef(mod)[2])

mm <- princomp(X.s[complete.cases(X.s),])
plot(mm)



# question 1
X <- cbind(X2, X3)
dd <- dist(X)
dd.h <- hclust(dd)
plot(dd.h)


# silhouette
par(mfrow=c(2,2))
for (i in 2:5)
 plot(silhouette(cutree(dd.h,i), dd))



# question 2
groups <- as.factor(cutree(dd.h, 3))
table(groups)
table(groups, X1$tipocancer)


XX <- X[groups!="3",]
dd <- dist(XX)
dd.h <- hclust(dd)
plot(dd.h, labels=FALSE)



# question 3
mod <- glm(casoc ~ groups + peso + sexo + Colesterol, 
           data=X1, family="binomial")
summary(mod)

# question 4
library(made4)
library(impute)
X.imp <- impute.knn(t(X))$data
mod <- ord(X.imp, classvec = groups)
plot(mod, nlab=3)

# Extra

library(compareGroups)
X$groups <- as.factor(groups)
mm <- compareGroups(groups ~ ., X)
createTable(mm)




