# Post-omic data analysis

4-hours course introducing post omic data analysis and visualization. 

**Lecturer**: Juan R Gonzalez, Associate Research Professor, Head of Bioinformatics Research Group in Genetic Epidemiology (BRGE) of
Barcelona Institute for Global Health (ISGlobal) and Adjunct Professor to the Mathematics 
Department of Universitat Autonoma de Barcelona (UAB)

- **url:** [BRGE](http://brge.isglobal.org)
- **GitHub:** [GitHub BRGE](https://github.com/isglobal-brge)
- **email:** juanr.gonzalez@isglobal.org

## Introduction

The ultimate goal of omic association studies is to find biomarkers that can be used in personalized medicine. Several methods (GWAS, EWAS, penalized regression, machine learning algorithms, ….) are used towards this end. However, it is essential to know whether the selected biomarkers (SNPs, transcripts, CpGs, …) have a functional impact, are involved in key genes or biological pathways. Post-omic data analysis performs an integrated visualization new experimental data and known genomic information. In this short course we will first illustrate how to visualize omic results to help data analysts and biologists to interpret the results obtained after association analyses. This will include how to create genomic circos plots which help, for instance, the visualization of interaction between variants, or how to visualize cis- and trans- transcriptome or epigenome effects. We will also illustrate how to create plot tracks using GViz to incorporate gene annotation and other features such as eQTLs, GWAS hits or motifs in Manhattan or Locus Zoom plots. Biological interpretation will be facilitated by describing how to annotate variants into databases from large genomic projects such as ENCODE, the 1000 Genomes, Roadmap Epigenomics. Enrichment analysis at pathway, disease, or environment level using different Bioconductor packages will also be shown. The methodology will consist on a short introduction describing the main methods followed by an illustration using real datasets. Then, the attendees will practice by analyzing real data problems. The R code, exercises and answers will be available at this GitHub repository.


## Outline

- **Session 1:**	
     - Introduction to Bioconductor
     - GenomicRanges
     - Bioconductor infrastructures (ExpressionSet and SummarizedExperiments, ...)

- **Session 2:**
     - Retreiving omic data annotation
     - Experimental data visualization and annotation using GViz 
     - Circos plots
     
- **Session 3:**
     - Enrichment analysis with curated databases (GO, KEGG, MSigDB, DISGeNET, ...)
     - Enrichment with public databases (ENCODE, ROADMAP Epigenomics, …)


## Material

- R code given in the slides are available at folder called `R`. 

- Data for reproducing the slides are available at folder called `data`.

- Data for the exercises are available at folder `data_exercises`.


## Install required packages

The material has been created using R version 3.5.0. These packages are required to be installed before starting the course:

```
install.packages("BiocManager")

library(BiocManager)
install(c("GenomicRanges", "GEOquery", "GenomicRanges",
          "SummarizedExperiment", "airway", "org.Hs.eg.db",
          "Gviz", "OmicCircos", "biomaRt", "GenomicFeatures",
          "TxDb.Hsapiens.UCSC.hg19.knownGene", "GO.db", 
          "AnnotationHub", "limma", "MEAL", "edgeR", 
          "GOstats", "KEGG.db", "clusterProfiler",
          "regioneR", "omicade4", "curatedTCGAData",
          "MultiAssayExperiment",
          "BSgenome.Hsapiens.UCSC.hg19.masked")
```