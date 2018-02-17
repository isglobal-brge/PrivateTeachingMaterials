#
# Task 1
#

pp <- "c:/juan/CREAL/GitHub/biomarkers_multiple_tables/"
load(file.path(pp, "data_exercises/diet.Rdata"))
X <- cbind(X1[, c("casoc", "estudios", "peso", "mets_10a", "Colesterol")],
           X3)
X.c <- X[complete.cases(X),]
mydat <- X.c[, c("Colesterol", "gra_procmeat", "gra_oilolives", 
                 "gra_vegetables", "gra_fruit", "casoc")]
mydists <- list(Colesterol="binomial",
                gra_procmeat="gaussian",
                gra_oilolives="gaussian",
                gra_vegetables="gaussian",
                gra_fruit="gaussian",
                casoc="binomial")

mydag <- ggm::DAG(Colesterol~casoc+gra_procmeat, gra_vegetable~casoc+gra_fruit,
                  gra_oilolives ~ casoc, gra_procmeat ~ casoc)
mydag

model <- fitabn(dag.m=mydag, data.df = mydat, data.dists = mydists, 
                 centre=TRUE, create.graph = TRUE)

model$mlik

plot(model$graph)


#
# Task 2
#

pp <- "c:/juan/CREAL/GitHub/biomarkers_multiple_tables/"
load(file.path(pp, "data_exercises/diet.Rdata"))

vars1 <- c("casoc", "sexo", "peso", "mets_10a", "Colesterol")
X <- cbind(X1[, vars1],  X3[,1:6])
mydat <- X[complete.cases(X),]
names(mydat)

dd <- c("binomial", "binomial", "gaussian", "gaussian", "binomial",
        rep("gaussian", 6))
dd
names(dd) <- c(vars1, names(X3)[1:6])
mydists <- as.list(dd)
mydists

ban <- retain <- matrix(0, nrow=ncol(mydat), ncol=ncol(mydat))
colnames(ban) <- rownames(ban) <- names(mydat)
colnames(retain) <- rownames(retain) <- names(mydat)

ii <- rep(5, ncol(mydat))
names(ii) <- names(mydat)
max.par <- as.list(ii)


mycache<-buildscorecache(data.df=mydat, data.dists=mydists,
                         dag.banned=ban, dag.retained=retain,
                         max.parents = max.par)
mp.dag<-mostprobable(score.cache=mycache)
mod.best <- fitabn(dag.m=mp.dag, data.df=mydat,
                   data.dists=mydists, create.graph=TRUE)

mod.best$mlik
graph::plot(mod.best$graph)
