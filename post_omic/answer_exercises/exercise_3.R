library(MEAL)
library(clusterProfiler)
library(org.Hs.eg.db)


load(file="data_exercises/GSE18123.Rdata")

ans.sva <- runPipeline(gse18123,
                       variable = "group",
                       sva = TRUE)

fit <- getProbeResults(ans.sva, coef=2, 
                       fNames=c("ID", "Gene Symbol"))

head(fit)

de <- fit[fit$adj.P.Val<0.05,]  # 5% FDR
de
deGenes <- de$`Gene Symbol`
deGenes

# .... gene symbol into Entrez (required)
deGenes <- deGenes[deGenes!=""]
deGenes <- unlist(mget(deGenes[1:5], envir=org.Hs.egSYMBOL2EG,
                       ifnotfound = NA))


# GO
ans.go <- enrichGO(gene = deGenes, ont = "BP",
                   OrgDb ="org.Hs.eg.db",
                   readable=TRUE,
                   pvalueCutoff = 0.05)
ans.go
as.data.frame(ans.go)


# KEGG
ans.kegg <- enrichKEGG(gene = deGenes,
                       organism = 'hsa',
                       pvalueCutoff = 0.05)
ans.kegg
as.data.frame(ans.kegg)


# DisGENET

gda <- read.delim("data/curated_gene_disease_associations.tsv.gz")
disease2gene <- gda[, c("diseaseId", "geneId")]
disease2name <- gda[, c("diseaseId", "diseaseName")]
ans.dis <- enricher(deGenes, TERM2GENE=disease2gene,
                    TERM2NAME=disease2name)
ans.dis
as.data.frame(ans.dis)
