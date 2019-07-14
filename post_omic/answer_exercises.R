
load(file="data_exercises/breast_tcga.RData")

roi <- GRanges(seqnames="chr6", IRanges(start=151.2e6,
                                 end=151.8e6))
breast.roi <- subsetByOverlaps(breast, roi)
assay(breast.roi)
rowRanges(breast.roi)


load(file="data_exercises/GSE85426.Rdata")


library(MEAL)
ans <- runPipeline(gse18123,
                   variable = "group")
ans
plot(ans, type="QQ")

ans.sva <- runPipeline(gse18123,
                   variable = "group",
                   sva = TRUE)
ans.sva

plot(ans.sva, coef=2, type="QQ")

# AUTISM
fit <- getProbeResults(ans.sva, coef=2, 
                       fNames=c("ID", "Gene Symbol"))

head(fit)


library(org.Hs.eg.db)
cols <- c("SYMBOL", "CHR", "CHRLOC", "CHRLOCEND"  )
annot <- select(org.Hs.eg.db, keys="KIAA1468",
                columns=cols, 
                keytype="SYMBOL")
annot

library(Gviz)
from <- annot$CHRLOC - 10e6
to <- annot$CHRLOC + 10e6
region <- GRanges(seqnames = "chr18", 
                  IRanges(start=from, end=to))
knownGenes <- UcscTrack(genome="hg19", chromosome="chr18",
                        track="knownGene", from=from, to=to,
                        trackType="GeneRegionTrack",
                        rstarts="exonStarts",
                        rends="exonEnds", gene="name",
                        symbol="name", transcript="name",
                        strand="strand", fill="#8282d2",
                        name="UCSC Genes")

cpgIslands <- UcscTrack(genome="hg19", chromosome="chr18",
                        track="cpgIslandExt", from=from,
                        to=to, trackType="AnnotationTrack",
                        start="chromStart", end="chromEnd",
                        id="name", shape="box",
                        fill="#006400", name="CpG Islands")

axTrack <- GenomeAxisTrack()
idxTrack <- IdeogramTrack(genome="hg19", chromosome="chr18")

plotTracks(list(idxTrack, axTrack, knownGenes, cpgIslands),
           from=from, to=to, showTitle=FALSE)


library(AnnotationHub)
ah <- AnnotationHub()

H3K27Ac <- query(ah , c("EpigenomeRoadMap", "H3K27Ac"))
H3K27Ac
H3K27Ac[grep("Cingulate_Gyrus", H3K27Ac$title),]
H3K27Ac[grep("Hippocampus", H3K27Ac$title),]
H3K27Ac[grep("Mid_Frontal_Lobe", H3K27Ac$title),]
H3K27Ac[grep("Adipose_Tissue", H3K27Ac$title),]


ids <- c("AH43527", "AH43543", "AH43571", "AH44156")
peaks <- list()
for (i in ids){
  peaks[[i]] <- H3K27Ac[[i]]
}

peaksTrack <- lapply(peaks,
                     DataTrack, data="score",
                     type="mountain", ylim=c(0,150),
                     name="H3K27Ac")



library(gwascat)
data(ebicat37UCSC) # genome tags and UCSC seqnames
gwas <- gwcex2gviz(basegr = ebicat37UCSC,
                       contextGR=region, plot.it=FALSE)

gwasTrack <- gwas[[1]]
genesTrack <- gwas[[3]]

plotTracks(c(idxTrack, axTrack, genesTrack, knownGenes, cpgIslands,
                peaksTrack, gwasTrack),
           from=from, to=to, showTitle=FALSE)


load(file="data_exercises/breast_tcga.RData")

library(edgeR)
counts <- assay(breast)
group <- breast$er

d <- DGEList(counts=airCounts, group=group)
design <- model.matrix( ~ group)
d <- estimateDisp(d, design)
d$common.dispersion
ans <- runPipeline(breast,
                   variable = "er",
                   sva = TRUE)
fit <- getProbeResults(ans, coef=4, 
                       fNames=c("ID", "Gene Symbol"))
head(fit)

pData(gse18123)$group <- factor(
pData(gse18123)$"diagnosis:ch1", 
levels = c("CONTROL", "AUTISM", "PDD-NOS", 
           "ASPERGER'S DISORDER"))
                                            )