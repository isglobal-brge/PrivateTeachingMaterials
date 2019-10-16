The complete material of this course (R code, datasets and html files with
the theory and examples) are available in this GitHub repository.

If you want to reproduce the examples and do the exercises, all existing
files should be downloaded into your lab top before starting the 
course. 

It is also highly recommended you install the following packages in your 
computer (R version 3.5.0 is enough). They are heavy packages that can 
take long time to be downloaded.


install.packages("BiocManager")

library(BiocManager)
install(c("GenomicRanges", "GEOquery", "GenomicRanges",
          "SummarizedExperiment", "airway", "org.Hs.eg.db",
          "Gviz", "OmicCircos", "biomaRt", "GenomicFeatures",
          "TxDb.Hsapiens.UCSC.hg19.knownGene", "GO.db", 
          "AnnotationHub", "limma", "MEAL", "edgeR", 
          "GOstats", "KEGG.db", "clusterProfiler",
          "regioneR", "gwascat", "rtracklayer",
          "BSgenome.Hsapiens.UCSC.hg19.masked")