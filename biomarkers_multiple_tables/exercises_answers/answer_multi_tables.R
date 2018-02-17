#
# Task 1
#

library(omicade4)
setwd("c:/juan/CREAL/GitHub/biomarkers_multiple_tables/")


load("data_exercises/diet.Rdata")
sel <- complete.cases(X2) & complete.cases(X3)
X2.comp <- scale(X2[sel, -which(names(X2)=="t_ethanol")])
X3.comp <- scale(X3[sel,])
X <- list(nut=t(X2.comp), food=t(X3.comp))
mm <- mcia(X)
plot(mm, df.col=c("red", "blue"), gene.nlab = 2)

library(made4)
XX <- cbind(X2, X3)[sel,]
XX <- scale(XX)
mm2 <- ord(t(XX))
plotgenes(mm2, nlab = 4)


#
# Task 2
#

# Question 1
setwd("c:/Juan/CREAL/GitHub/biomarkers_multiple_tables/")
load("data_exercises/nci60.Rdata")

# Question 2 
library(omicade4)
ans.mcia <- mcia(nci60, cia.nf = 2)

plot.mcia(ans.mcia, sample.lab = TRUE,  
          phenovec = cancer, 
          gene.nlab = 5, sample.legend = FALSE,
          df.color =c("cyan", "magenta", "red4"))

# these are associated with Leukemia
topVar(ans.mcia, axis=1, topN=5, end="pos")

# these are associated with Melanoma
topVar(ans.mcia, axis=1, topN=5, end="neg")

# these are associated with Central Nervous System
topVar(ans.mcia, axis=2, topN=5, end="neg")


# these are associated with Melanoma
topVar(ans.mcia, axis=2, topN=5, end="pos")



#
# TASK 3
#

library(RGCCA)
load("data_exercises/nci60.Rdata")
source("Day3-integration_multiple_tables/R/plotInd.R")
source("Day3-integration_multiple_tables/R/selectVars.R")



X1 <- t(nci60$mrna)
X2 <- t(nci60$miRNA)
X3 <- t(nci60$prot)
Y <- model.matrix( ~ cancer)[,-1]
A <- list(X1, X2, X3, Y)
lapply(A, dim)

C <-  matrix(c(0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0), nrow=4, ncol=4)
C

sgcca.nci60 = sgcca(A, C, c1 = c(.071, .3, 0.1, 1),
                    ncomp = c(2, 2, 2, 2),
                    scheme = "centroid",
                    scale = TRUE,
                    verbose = FALSE)

plotInd(sgcca.nci60, cancer)

# miRNA Leukemia
selectVars(sgcca.nci60, table = 2, axis = 1, end = "pos")
# miRNA Melanoma
selectVars(sgcca.nci60, table = 2, axis = 1, end = "neg")
# miRNA CNS
selectVars(sgcca.nci60, table = 2, axis = 2, end = "neg")


#
# Shrinkage estimation
#

p1 <- tau.estimate(X1)
p2 <- tau.estimate(X2)
p3 <- tau.estimate(X3)

sgcca.nci60 <- sgcca(A, C, c1 = c(p1, p2, p3, 1),
                    ncomp = c(2, 2, 2, 2),
                    scheme = "centroid",
                    scale = TRUE,
                    verbose = FALSE)

plotInd(sgcca.nci60, cancer)
length(selectVars(sgcca.nci60, table=1, axis=1, end="pos"))
length(selectVars(sgcca.nci60, table=1, axis=1, end="neg"))
length(selectVars(sgcca.nci60, table=1, axis=2, end="neg"))

sgcca.nci60 <- sgcca(A, C, c1 = c(0.06, p2, p3, 1),
                     ncomp = c(2, 2, 2, 2),
                     scheme = "centroid",
                     scale = TRUE,
                     verbose = FALSE)



###
###  Nutrients
###

load("data_exercises/diet.Rdata")
sel <- complete.cases(X2) & complete.cases(X3) & complete.cases(X1$casoc)
X2.comp <- scale(X2[sel, -which(names(X2)=="t_ethanol")])
X3.comp <- scale(X3[sel,])
Y <- model.matrix(~X1$casoc[sel])[,-1]
X <- list(nut=t(X2.comp), food=t(X3.comp), casoc=Y)


C <-  matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), nrow=3, ncol=3)
C


sgcca.nci60 <- sgcca(X, C, c1 = c(0.6, 0.6, 1),
                     ncomp = c(1, 1, 1),
                     scheme = "centroid",
                     scale = TRUE,
                     verbose = FALSE)
Y
