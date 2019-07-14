#
# EXERCISE 2
#

load(file="data_exercises/GSE18123.Rdata")
table(gse18123$group)

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


de <- fit[fit$adj.P.Val<0.05,]  # 5% FDR
dim(de)

de$sign <- ifelse(de$logFC < 0, "red", "blue")
table(de$sign)

library(org.Hs.eg.db)
# see columns(org.Hs.eg.db)
cols <- c("SYMBOL", "CHR", "CHRLOC", "CHRLOCEND"  )
annot <- select(org.Hs.eg.db, keys=de$"Gene Symbol",
                columns=cols, 
                keytype="SYMBOL")
annot <- subset(annot, !duplicated(SYMBOL) & !is.na(CHR))
bd <- merge(annot, de, by.x="SYMBOL", by.y="Gene Symbol")
dim(bd)

library(OmicCircos)
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



library(org.Hs.eg.db)
cols <- c("SYMBOL", "CHR", "CHRLOC", "CHRLOCEND"  )
annot <- select(org.Hs.eg.db, keys="KIAA1468",
                columns=cols, keytype="SYMBOL")
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

