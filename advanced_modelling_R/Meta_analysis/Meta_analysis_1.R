## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, comment="", message=FALSE, warning=FALSE, cache=TRUE, fig.width = 4, fig.height = 4)
options(width=80)


## ----install, eval=FALSE-------------------------------------------------
## install.packages(c("rmeta", "meta", "altmeta",
##                    "RISmed", "RCurl"))


## ----load_library--------------------------------------------------------
library("RISmed")


## ------------------------------------------------------------------------
rofecoxib <- EUtilsSummary("rofecoxib[ti] +
                           British Medical Journal[jo]", 
                           db = "pubmed")


## ----exact---------------------------------------------------------------
QueryTranslation(rofecoxib)


## ----n-------------------------------------------------------------------
QueryCount(rofecoxib)


## ----metadata------------------------------------------------------------
metadata <- EUtilsGet(rofecoxib)
metadata # Medline Object


## ----methods-------------------------------------------------------------
ls("package:RISmed")


## ----example1------------------------------------------------------------
ArticleTitle(metadata)[1:5]


## ----example2------------------------------------------------------------
Author(metadata)[[1:2]]


## ----example3------------------------------------------------------------
YearPubmed(metadata)


## ----plotYear0, eval=FALSE-----------------------------------------------
## pm <- EUtilsSummary("viagra[ti]",
##                     db = "pubmed")
## metadata.pm <- EUtilsGet(pm)
## y <- YearPubmed(metadata.pm)
## hist(y, ylab = "Number of articles",
##      xlab = paste0("Query date: ",
##                   Sys.Date()),
##      main = paste0("PubMed articles containing
##                    viagra = ", length(y)))


## ----plotYear, echo=FALSE------------------------------------------------
pm <- EUtilsSummary("viagra[ti]",
                    db = "pubmed")
metadata.pm <- EUtilsGet(pm)
y <- YearPubmed(metadata.pm)
hist(y, ylab = "Number of articles",
     xlab = paste0("Query date: ", 
                  Sys.Date()),
     main = paste0("PubMed articles containing 
                   viagra = ", length(y)))


## ----advance-------------------------------------------------------------
AuthorList <- Author(metadata) # Extract list of authors
LastFirst <- sapply(AuthorList, function(x) paste(x$LastName,
                                                  x$ForeName))
sort(table(unlist(LastFirst)), dec = TRUE)[1:3] # Tabulate & Sort


## ----create_file, echo=FALSE---------------------------------------------
data("woodyplants", package="meta")
study <- paste("Study", 1:10)
temp <- subset(woodyplants, treat=="fert")[1:10, ]
example <- data.frame(Study = study,
  Ne=temp$n.elev, Me=temp$mean.amb, Se=temp$sd.amb,
  Nc=temp$n.elev, Mc=temp$mean.elev, Sc=temp$sd.elev)
write.table(example, file="../data/example.txt",
            quote=FALSE, sep="\t", row.names = FALSE)
example


## ----amlopidine----------------------------------------------------------
library(meta)
data("amlodipine", package="meta")
res.cont <- metacont(n.amlo, mean.amlo, sqrt(var.amlo),
                     n.plac, mean.plac, sqrt(var.plac),
                     data=amlodipine, studlab=study)


## ----load_data-----------------------------------------------------------
data("cochrane", package="rmeta")
cochrane


## ----meta_cochran--------------------------------------------------------
library(rmeta)
cochrane
model.FE <- meta.MH(n.trt, n.ctrl, ev.trt, ev.ctrl, 
                    names=name, data=cochrane)
model.RE <- meta.DSL(n.trt, n.ctrl, ev.trt, ev.ctrl,
                     names=name, data=cochrane)


## ----meta_cochran_res----------------------------------------------------
model.FE 
model.RE 


## ----anal_rma------------------------------------------------------------
library(meta)
res.bin.FE <- metabin(ev.trt, n.trt, ev.ctrl, n.ctrl, 
               data=cochrane, studlab = name, 
               comb.random = FALSE)
res.bin.RE <- metabin(ev.trt, n.trt, ev.ctrl, n.ctrl, 
               data=cochrane, studlab = name, 
               comb.fixed = FALSE)


## ----sum_rma-------------------------------------------------------------
summary(res.bin.FE)
summary(res.bin.RE)


## ----sum_rma_both--------------------------------------------------------
res.bin <- metabin(ev.trt, n.trt, ev.ctrl, n.ctrl, 
                   data=cochrane, studlab = name)
res.bin


## ----contrib-------------------------------------------------------------
vari <- (res.bin$seTE)^2
contrib <- 1/vari/sum(1/vari) * 100
barplot(contrib, names = res.bin$studlab, ylim = c(0, 60), 
        las = 2, col = "royalblue")


## ----diffe_var-----------------------------------------------------------
estimators <- c("DL", "REML", "HE", "HS", "SJ", "ML", "EB")
taus <- sapply(estimators, function(method) {
  metabin(ev.trt, n.trt, ev.ctrl, n.ctrl, 
          data=cochrane, method.tau=method)$tau
})


## ----plot_diffe_var,eval=FALSE-------------------------------------------
## barplot(taus, las=2, col="tomato",
##         ylab=expression(tau^2))


## ----plot_diffe_var2, echo=FALSE-----------------------------------------
barplot(taus, las=2, col="tomato", 
        ylab=expression(tau^2))


## ----get_q---------------------------------------------------------------
res.bin$Q
res.bin$pval.Q


## ----example-------------------------------------------------------------
res.bin


## ----plot_rma, fig.width=10----------------------------------------------
forest(res.bin, studlab=cochrane$name)


## ----plot2_rma,  fig.width=10--------------------------------------------
forest(res.bin, layout = "JAMA")


## ----plot3_rma, fig.width=10---------------------------------------------
forest(res.bin, sortvar=-TE, comb.fixed=FALSE)


## ----prediction, fig.width=10--------------------------------------------
forest(res.bin, prediction = TRUE)


## ----getBetas------------------------------------------------------------
studies <- paste0("Study", 1:5)
or <- c(1.23, 1.36, 1.08, 1.24, 1.81)
ciInf <- c(1.12, 1.06, 0.97, 1.20, 1.55)
beta <- log(or)
beta
se <- (beta - log(ciInf)) / 1.96 
se


## ----plotBetas, fig.width=10---------------------------------------------
res.obs <- metagen(beta, se)
forest(res.obs, layout = "JAMA")


## ----plotBetas2, fig.width=10--------------------------------------------
res.obs <- metagen(beta, se, sm="OR")
forest(res.obs, layout = "JAMA")


## ----sensit--------------------------------------------------------------
lev1out <- metainf(res.bin, pooled = "random")
lev1out


## ----altmeta-------------------------------------------------------------
library(altmeta)
out <- metaoutliers(res.cont$TE,      # observed effect sizes
                    res.cont$seTE^2)  # within-studies variances
out$outliers


## ------------------------------------------------------------------------
sessionInfo()

