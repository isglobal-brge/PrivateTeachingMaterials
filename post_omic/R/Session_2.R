## ----style, echo = FALSE, results = 'asis'-------------------------------
options(width=120)
knitr::opts_chunk$set(cache=TRUE, fig.align = TRUE, 
                      warning = FALSE,
                      message = FALSE, comment = "")


## ----ucsc_genes----------------------------------------------------------
library(Gviz)
from <- 65921878
to <- 65980988
knownGenes <- UcscTrack(genome="mm9", chromosome="chrX",
                        track="knownGene", from=from, to=to,
                        trackType="GeneRegionTrack",
                        rstarts="exonStarts",
                        rends="exonEnds", gene="name",
                        symbol="name", transcript="name",
                        strand="strand", fill="#8282d2",
                        name="UCSC Genes")


## ----other_refs----------------------------------------------------------
library(Gviz)
refGenes <- UcscTrack(genome="mm9", chromosome="chrX",
                      track="xenoRefGene", from=from, to=to,
                      trackType="GeneRegionTrack",
                      rstarts="exonStarts", rends="exonEnds",
                      gene="name",
                      symbol="name2", transcript="name",
                      strand="strand", fill="#8282d2",
                      stacking="dense", name="Other RefSeq")

ensGenes <- UcscTrack(genome="mm9", chromosome="chrX",
                      track="ensGene", from=from, to=to,
                      trackType="GeneRegionTrack",
                      rstarts="exonStarts", rends="exonEnds",
                      gene="name",
                      symbol="name2", transcript="name",
                      strand="strand", fill="#960000",
                      name="Ensembl Genes")


## ----annot_CpG_SNP-------------------------------------------------------
cpgIslands <- UcscTrack(genome="mm9", chromosome="chrX",
                        track="cpgIslandExt", from=from,
                        to=to, trackType="AnnotationTrack",
                        start="chromStart", end="chromEnd",
                        id="name", shape="box",
                        fill="#006400", name="CpG Islands")

snpLocations <-  UcscTrack(genome="mm9", chromosome="chrX",
                           track="snp128", from=from, to=to,
                           trackType="AnnotationTrack",
                           start="chromStart", end="chromEnd",
                           id="name",  feature="func",
                           strand="strand", shape="box",
                           stacking="dense", fill="black",
                           name="SNPs")


## ----annot_conservation_CpG----------------------------------------------
conservation <- UcscTrack(genome="mm9", chromosome="chrX",
                          track="Conservation",
                          table="phyloP30wayPlacental",
                          from=from, to=to,
                          trackType="DataTrack",
                          start="start", end="end",
                          data="score",
                          type="hist", window="auto",
                          col.histogram="darkblue",
                          fill.histogram="darkblue",
                          ylim=c(-3.7, 4),
                          name="Conservation")

gcContent <- UcscTrack(genome="mm9", chromosome="chrX",
                       track="GC Percent", table="gc5Base",
                       from=from, to=to,
                       trackType="DataTrack", start="start",
                       end="end", data="score",
                       type="hist", window=-1,
                       windowSize=1500,
                       fill.histogram="black",
                       col.histogram="black",
                       ylim=c(30, 70), name="GC Percent")


## ----tables_UCSC---------------------------------------------------------
library(rtracklayer) 
session <- browserSession() 
genome(session) <- "hg19" 
trackNames(session)[1:20]


## ----genome_ideogram_chrX------------------------------------------------
axTrack <- GenomeAxisTrack()
idxTrack <- IdeogramTrack(genome="mm9", chromosome="chrX")


## ----plot_UCSC-----------------------------------------------------------
plotTracks(list(idxTrack, axTrack, knownGenes, refGenes,
                ensGenes, cpgIslands,
                gcContent, conservation, snpLocations),
           from=from, to=to, showTitle=FALSE)


## ----gwascat-------------------------------------------------------------
library(gwascat)
data(ebicat37) # gwas catalog - h19 coordinates
roi <- GRanges(seqnames="chr8", 
               IRanges(start = 48e6, end = 56e6))

gwas <- gwcex2gviz(basegr = ebicat37, 
                   contextGR=roi, plot.it=FALSE)
gwasTrack <- gwas[[1]]
gwasTrack
genesTrack <- gwas[[3]]
genesTrack


## ----plot_example_GWAS---------------------------------------------------
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = "hg19", 
                        chromosome = "chr8")

plotTracks(list(itrack, gtrack, genesTrack,
                gwasTrack), 
           from = start(roi),
           to = end(roi))


## ----load_h3k4me3, eval=FALSE--------------------------------------------
## H3K4me3 <- UcscTrack(genome = "hg19",
##                      chromosome = seqnames(region),
##                      track = "Layered H3K27Ac",
##                      table = "wgEncodeBroadHistoneHsmmH3k4me3StdSig",
##                      from=1, to=10e6,
##                      trackType = "DataTrack",
##                      start = "start",
##                      end = "end", data="score", type="hist",
##                      window=-1, windowSize=1500,
##                      fill.histogram="orange",
##                      col.histogram="orange", name="H3K27Ac")


## ----roadmap-------------------------------------------------------------
library(AnnotationHub)
ah <- AnnotationHub()

H3K27Ac <- query(ah , c("EpigenomeRoadMap", "H3K27Ac"))
H3K27Ac


## ----look_ids_roadmap----------------------------------------------------
H3K27Ac[grep("Lung", H3K27Ac$title),]
H3K27Ac[grep("Stomach", H3K27Ac$title),]
H3K27Ac[grep("Bladder", H3K27Ac$title),]


## ----get_ke--------------------------------------------------------------
ids <- c("AH44566", "AH44076", "AH44180")
peaks <- list()
for (i in ids){
 peaks[[i]] <- H3K27Ac[[i]]
}

peaksTrack <- lapply(peaks, DataTrack, data="score",
                     type="mountain", ylim=c(0,100),
                     name="H3K27Ac")


## ----load_pleio----------------------------------------------------------
library(GenomicRanges)
rois <- read.delim("data/ng.3570-S2.txt", comment.char = "#",
                  as.is=TRUE)
colnames(rois)
regionsGR <- makeGRangesFromDataFrame(rois, seqnames="chr", 
                                      start.field = "st",
                                      end.field = "sp")
regionsGR


## ----roi_pleio_chr1------------------------------------------------------
rr <- GRanges(seqnames = "chr1", 
              IRanges(start = 1e6, end=100e6))
rr2 <- subsetByOverlaps(regionsGR, rr)
roi <- AnnotationTrack(rr2, name="ROI")


## ----plot_pleio----------------------------------------------------------
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = "hg19", 
                        cromosome=seqnames(rr2))


plotTracks(c(itrack, gtrack, roi, peaksTrack), 
           from = min(start(rr2)),
           to = max(end(rr2)))


## ----add_genes_pleio-----------------------------------------------------
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
allg <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
allg.range <- subsetByOverlaps(allg, rr2)
allg.range$symbol <- mapIds(Homo.sapiens::Homo.sapiens, 
                                keys=allg.range$gene_id,
                                keytype="ENTREZID",
                                column="SYMBOL")

genes <- GeneRegionTrack(allg.range, genome = "hg19",
                         chromosome = seqnames(rr2),
                         showId=TRUE,
                         geneSymbol=TRUE,
                         start = min(start(rr2)), 
                         end = max(end(rr2)),
                         name = "Genes")


## ----plot_pleio_genes, fig.height=14, fig.width=8------------------------
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = "hg19", 
                        cromosome=seqnames(rr2))

plotTracks(c(itrack, gtrack, roi, genes, peaksTrack), 
           from = min(start(rr2)),
           to = max(end(rr2)))


## ----de------------------------------------------------------------------
library(limma)
load("data/GSE69683.Rdata")
design <- model.matrix( ~ characteristics_ch1,
                        data=gse69683.expr)
colnames(design)
fit <- lmFit(gse69683.expr, design)
fit <- eBayes(fit)
de <- topTable(fit, coef=4)
de[1:5 , c(17:22)]


## ----qqplot_asthma-------------------------------------------------------
qqt(fit$t, df=fit$df.prior+fit$df.residual,
    pch=16,cex=0.2)
abline(0,1, col="red", lwd=2)


## ----de_meal-------------------------------------------------------------
library(MEAL)
load("data/GSE69683.Rdata")

# shows available annotated variabels
names(fData(gse69683.expr))

# DE analysis
ans <- runPipeline(gse69683.expr,
                   variable = "characteristics_ch1",
                   sva = TRUE)
ans

fit <- getProbeResults(ans, coef=4, 
                       fNames=c("ID", "Gene Symbol"))
head(fit)


## ----qqplot_meal---------------------------------------------------------
plot(ans, type="QQ")


## ----de_10fdr------------------------------------------------------------
de <- fit[fit$adj.P.Val<0.1,]
dim(de)
de


## ----down_up-------------------------------------------------------------
de$sign <- ifelse(de$logFC < 0, "red", "blue")
table(de$sign)


## ----circos--------------------------------------------------------------
library(org.Hs.eg.db)
# see columns(org.Hs.eg.db)
cols <- c("SYMBOL", "CHR", "CHRLOC", "CHRLOCEND"  )
annot <- select(org.Hs.eg.db, keys=de$"Gene Symbol",
                columns=cols, 
                keytype="SYMBOL")
annot <- subset(annot, !duplicated(SYMBOL) & !is.na(CHR))
bd <- merge(annot, de, by.x="SYMBOL", by.y="Gene Symbol")
dim(bd)
head(bd)


## ----omic1, fig.height=5, fig.width=5------------------------------------
library(OmicCircos)

plot(c(1,800), c(1,800), 
     type="n", axes=FALSE, xlab="", ylab="")
circos(R=300, type="chr", cir="hg19", col=TRUE, 
       print.chr.lab=TRUE, W=10, cex=2)


## ----omic2, fig.height=5, fig.width=5------------------------------------
plot(c(1,800), c(1,800), 
     type="n", axes=FALSE, xlab="", ylab="")
circos(R=300, type="chr", cir="hg19", col=TRUE, 
       print.chr.lab=TRUE, W=10, cex=2)

circos(R=250, cir="hg19", W=40,
      mapping = bd[,c("CHR", "CHRLOCEND", "logFC")], 
      type="b3", lwd=3, B=TRUE, cutoff = 0, 
      col=bd$sign)

circos(R=300, cir="hg19", W=40, 
       mapping=bd[,c("CHR", "CHRLOCEND", "SYMBOL")], 
       type="label", side="out", cex=0.6)


## ----table_grouping_var--------------------------------------------------
table(gse69683.expr$`characteristics_ch1`)


## ----de_meal_2_groups----------------------------------------------------
ff <- list()
for (i in 1:2)
  ff[[i]] <- getProbeResults(ans, coef=i+2, 
                       fNames=c("ID", "Gene Symbol"))
de <- lapply(ff, function(x) x[x$adj.P.Val<0.1,])
de <- lapply(de, function(x) cbind(x, 
                                   sign=ifelse(x$logFC < 0,
                                               "red",
                                                "blue")))
de <- lapply(de, function(x) merge(annot, x, 
                                   by.x="SYMBOL", 
                                   by.y="Gene Symbol"))
head(de[[1]])


## ----circos_three, fig.height=5, fig.width=5-----------------------------

plot(c(1,800), c(1,800), type="n", 
     axes=FALSE, xlab="", ylab="")

circos(R=300, type="chr", cir="hg19", col=TRUE, 
       print.chr.lab=TRUE, W=10, cex=2)

mycols <- c("green", "darkblue")
for (i in 1:2) {
  circos(R=300, cir="hg19", W=40, 
       mapping=de[[i]][,c("CHR", "CHRLOCEND", "SYMBOL")], 
       type="label", side="out", cex=0.6)
  
  circos(R=250 - 40*(i-1), cir="hg19", W=40,
       mapping = de[[i]][,c("CHR", "CHRLOCEND", "logFC")], 
       type="b3", lwd=3, B=TRUE, cutoff = 0, 
       col=de[[i]]$sign)
}


## ------------------------------------------------------------------------
sessionInfo()

