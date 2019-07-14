## ----style, echo = FALSE, results = 'asis'-------------------------------
options(width=120)
knitr::opts_chunk$set(cache=TRUE, fig.align = TRUE, 
                      warning = FALSE,
                      message = FALSE, comment = "")


## ----install, eval=FALSE-------------------------------------------------
## install.packages("BiocManager")
## library(BiocManager)
## install(c("DESeq2", "org.Hs.eg.db"))
## 
## # or
## 
## BiocManager::install("DESeq2")


## ----install_github, eval=FALSE------------------------------------------
## install.packages("devtools")
## devtools::install_github("isglobal-brge/SNPassoc")


## ----require-------------------------------------------------------------
library(GenomicRanges)


## ----help-bioc, eval=FALSE-----------------------------------------------
## help(package="GenomicRanges")
## vignette(package="GenomicRanges")
## vignette(package="GenomicRanges",
##          "GenomicRangesHOWTOs")
## ?GRanges


## ----download_geo, eval=FALSE--------------------------------------------
## library(GEOquery)
## gse69683 <- getGEO("GSE69683", destdir = ".")
## gse69683.expr <- gse69683[[1]]


## ----load_geo------------------------------------------------------------
load("data/GSE69683.Rdata")


## ----show_geo------------------------------------------------------------
gse69683.expr


## ----exprs---------------------------------------------------------------
expr <- exprs(gse69683.expr)
dim(expr)
expr[1:5,1:5]


## ----pheno---------------------------------------------------------------
pheno <- phenoData(gse69683.expr)
pheno
colnames(pheno)[1:10]


## ----show----------------------------------------------------------------
group <- pheno$characteristics_ch1
table(group)


## ----boxplot-------------------------------------------------------------
boxplot(expr["1007_PM_s_at",] ~ group)


## ----get_annot-----------------------------------------------------------
probes <- fData(gse69683.expr)
probes[1:5, 1:5]


## ----subset_gse----------------------------------------------------------
sel <- "cohort: Healthy, non-smoking"
mask <- gse69683.expr$characteristics_ch1%in%sel
gse <- gse69683.expr[ , mask]
gse


## ----createGR------------------------------------------------------------
library(GenomicRanges)
gr <- GRanges(seqnames=c(rep("chr1", 4), rep("chr2", 4)),
              ranges = IRanges(start = c(1000, 1800, 5300, 7900,
                                         1300, 2100, 3400, 6700),
                               end =c(2200, 3900, 5400, 8100,
                                      2600, 3300, 4460, 6850)),
              strand = rep(c("+", "-"), 4),
              disease = c(rep("Asthma",4), rep("Obesity",4)))
gr


## ----gr1-----------------------------------------------------------------
gr[1]


## ----gr2-----------------------------------------------------------------
seqnames(gr)
seqnames(gr)[1] <- "chr2"
gr


## ----gr3-----------------------------------------------------------------
gr$gene_id <- paste0("Gene", 1:8)
gr


## ----column_granges------------------------------------------------------
mcols(gr)
table(mcols(gr)$disease)


## ----gr4-----------------------------------------------------------------
#shift: move all intervals 10 base pair towards the end
shift(gr, 10)

#shift: move each intervals individually
shift(gr, seq(10,100, length=8))

#flank:  recover regions next to the input set. 
#        For a 50 base stretch upstream (negative value for
#        downstream)
flank(gr, 50)


## ----gr8-----------------------------------------------------------------
target <- GRanges(seqnames="chr1", 
                  range=IRanges(start=1200, 4000))
target
gr.ov <- findOverlaps(target, gr)
gr.ov


## ----airway--------------------------------------------------------------
library(SummarizedExperiment)
data(airway, package="airway")
se <- airway
se


## ----exp_data------------------------------------------------------------
names(assays(se))
gene.dat <- assays(se)$counts
gene.dat[1:5, 1:5]


## ----pheno_data----------------------------------------------------------
colData(se)


## ----pheno_data_treated--------------------------------------------------
se[, se$dex == "trt"]


## ----subset_interval-----------------------------------------------------
roi <- GRanges(seqnames="chr1", ranges=100000:1100000)
# or
roi <- GRanges(seqnames="1", IRanges(start=100000,
                                        end=1100000))
se.roi <- subsetByOverlaps(se, roi)
se.roi


## ----chr_style-----------------------------------------------------------
seqlevelsStyle(se)


## ------------------------------------------------------------------------
sessionInfo()

