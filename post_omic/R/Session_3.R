## ----style, echo = FALSE, results = 'asis'-------------------------------
options(width=120)
knitr::opts_chunk$set(cache=TRUE, fig.align = TRUE, 
                      warning = FALSE,
                      message = FALSE, comment = "")


## ----deTbl---------------------------------------------------------------
deTable <-  matrix(c(28, 142, 501, 12000),
            nrow = 2,
            dimnames = list(DE=c("yes","no"),
                            GeneSet=c("in","out")))
              
deTable


## ----fisher--------------------------------------------------------------
fisher.test(deTable, alternative = "greater")


## ----loadEBrowser--------------------------------------------------------
library(EnrichmentBrowser)


## ----pdataAirway---------------------------------------------------------
library(airway)
data(airway)
table(airway$dex)


## ----DE_airway-----------------------------------------------------------
library(edgeR)
airCounts <- assay(airway)
group <- airway$dex
d <- DGEList(counts=airCounts, group=group)
design <- model.matrix( ~ group)
d <- estimateDisp(d, design)
d$common.dispersion

# LRT test
fit <- glmFit(d, design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)


## ----de_search-----------------------------------------------------------
tt <- topTags(lrt, n=Inf)
mask <- tt$table$FDR < 0.01 &
        abs(tt$table$logFC) > log2(2)
deGenes <- rownames(tt$table[mask, ])
head(deGenes)
length(deGenes)


## ----gene_universe-------------------------------------------------------
geneUniverse <- rownames(tt$table)
length(geneUniverse)


## ----genes_to_entrez-----------------------------------------------------
library(org.Hs.eg.db)
deGenes <- unlist(mget(deGenes, envir=org.Hs.egENSEMBL2EG,
                       ifnotfound = NA))

geneUniverse <- unlist(mget(geneUniverse, envir=org.Hs.egENSEMBL2EG,
                       ifnotfound = NA))


## ----functional_air------------------------------------------------------
library(GOstats)
params <- new("GOHyperGParams", geneIds=deGenes,
              universeGeneIds=geneUniverse,
              annotation="org.Hs.eg.db", ontology="BP",
              pvalueCutoff=0.05, conditional=FALSE,
              testDirection="over")

hgOver <- hyperGTest(params)
hgOver


## ----summary_hyperG_test_over--------------------------------------------
head(summary(hgOver))


## ----report_hyperG_test_over, eval=FALSE---------------------------------
## htmlReport(hgOver, file="goBPuncond.html")


## ----kegg_enrich---------------------------------------------------------
library(KEGG.db)
params.kegg <- new("KEGGHyperGParams", geneIds=deGenes,
                   universeGeneIds=geneUniverse,
                   annotation="org.Hs.eg.db", 
                   pvalueCutoff=0.05, 
                   testDirection="over")

hgOver.kegg <- hyperGTest(params.kegg)
head(summary(hgOver.kegg))


## ----GO_clusterProfiler--------------------------------------------------
library(clusterProfiler)
ans.go <- enrichGO(gene = deGenes, ont = "BP",
                   OrgDb ="org.Hs.eg.db",
                   readable=TRUE,
                   pvalueCutoff = 0.05)
tab.go <- as.data.frame(ans.go)
tab.go<- subset(tab.go, Count>5)
tab.go[1:5, 1:6]


## ----KEGG_clusterProfiler------------------------------------------------
ans.kegg <- enrichKEGG(gene = deGenes,
                       organism = 'hsa',
                       pvalueCutoff = 0.05)
tab.kegg <- as.data.frame(ans.kegg)
tab.kegg<- subset(tab.kegg, Count>5)
tab.kegg[1:5, 1:6]


## ----disgenet_import-----------------------------------------------------
gda <- read.delim("data/curated_gene_disease_associations.tsv.gz")
disease2gene <- gda[, c("diseaseId", "geneId")]
disease2name <- gda[, c("diseaseId", "diseaseName")]


## ----disgenet_clusterProfiler--------------------------------------------
ans.dis <- enricher(deGenes, TERM2GENE=disease2gene,
                    TERM2NAME=disease2name)
tab.dis <- as.data.frame(ans.dis)
tab.dis<- subset(tab.dis, Count>5)
tab.dis[,1:6]


## ----c3_tf---------------------------------------------------------------
c3.tf <- read.gmt("c:/Juan/CREAL/HELIX/pathways/GSEA/c3.tft.v6.2.entrez.gmt")

ans.tf <- enricher(deGenes, TERM2GENE=c3.tf)
tab.tf <- as.data.frame(ans.tf)
tab.tf<- subset(tab.tf, Count>5)
tab.tf[1:5,1:5]


## ----oad_brca_sig--------------------------------------------------------
load("data/brcaSigCNV.Rdata")
brca.gr.sig


## ----load_hub, eval=FALSE------------------------------------------------
## library(AnnotationHub)
## ah <- AnnotationHub()
## ahDb <- query(ah, pattern = c("Homo sapiens", "EnsDb"))
## ahDb


## ----get_coord, eval=FALSE-----------------------------------------------
## ahEdb <- ahDb[["AH60977"]]
## hg.genes <- genes(ahEdb)


## ----load_hg_genes, echo=FALSE-------------------------------------------
load("data/hgGenes.Rdata")


## ----show_seq_names------------------------------------------------------
seqnames(hg.genes)
seqnames(brca.gr.sig)


## ----get_genes, eval=FALSE-----------------------------------------------
## seqlevels(hg.genes) <- paste0("chr", seqlevels(hg.genes))


## ----sel_chr1------------------------------------------------------------
sel.protein <- hg.genes[seqnames(hg.genes)%in%c("chr1")]
sel.protein <- sel.protein[sel.protein$gene_biotype ==
                             "protein_coding"]
sel.cnvs <- brca.gr.sig[seqnames(brca.gr.sig)%in%c("chr1")]
length(sel.cnvs)


## ----test_cnvs_enrich----------------------------------------------------
library(regioneR)
library(BSgenome.Hsapiens.UCSC.hg19.masked)
res.prot <- overlapPermTest(A=sel.cnvs, B=sel.protein, 
                            mask=NA, genome="hg19",
                            ntimes=100,
                            per.chromosome=TRUE,
                            count.once=TRUE)
res.prot


## ----enrich_promoters_get, eval=FALSE------------------------------------
## file.prom <- "http://gattaca.imppc.org/regioner/data/UCSC.promoters.hg19.bed"
## promoters <- toGRanges(file.prom)
## sel.prom <- promoters[seqnames(promoters) %in% c("chr1")]


## ----enrich_promoters_get2, echo=FALSE-----------------------------------
load("data/selProm.Rdata")


## ----enrich_promoters----------------------------------------------------
res.prom <- overlapPermTest(sel.cnvs, sel.prom, 
                            ntimes=100, genome="hg19", 
                            count.once=TRUE)
res.prom


## ----sessionInfo---------------------------------------------------------
sessionInfo()

