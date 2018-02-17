# Methods to integrate multiple tables in biomedical studies to detect biomarkers and stratify individuals

5-days course introducing conditional and multivariate methods to analyze one or multiple tables in disease association studies. 
Course organized at Centro Nacional de Epidemiologia, 25-29 September, 2017.


**Lecturer**: Juan R Gonzalez, Associate Research Professor, Head of Bioinformatics Research Group in Genetic Epidemiology (BRGE) of
Barcelona Institute for Global Health (ISGlobal) and Adjunct Professor to the Mathematics 
Department of Universitat Autonoma de Barcelona (UAB)

- **url:** [BRGE](http://brge.isglobal.org)
- **GitHub:** [GitHub BRGE](https://github.com/isglobal-brge)
- **email:** juanr.gonzalez@isglobal.org

## Introduction

Omic data analysis is receiveing lot ot attention not only in genomic or genetic studies, but also in other biomedical settings. For instance, epidemiological studies are generating different types of omic data appart from SNP arrays such as methylation or gene expression. This makes extremely relevant to know existing methodologies to properly analyze such type of information to address most of the current questions in health association studies. This course aims to introduce (in a very general way) main methodologies and bioinformatic tools to manage, integrate and analyze different datasets. The methods include multivariate techniques such as cannonical correlation, coinertia analysis, O2PLS and reguralized generalized canonical correlation. Conditional and mendelian randomization will also be presented.

## Methodology

- Lectures will introduce statistical methods and how to analyze real data (from different settigs including nutrition, air pollution, genetic, genomics and other biomedical datasets using R packages

- Excersises will analyze data from a case/control study in cancer setting were controls and 4 types of cancer (colorectal, stomach, breast and prostate) were studied
  
- Datasets will include a set of variables that are considered as confounding variables, another set about nutrients and a third one ecoding food compsumption variables

## Outline

- **Day 1:** Introduction to R
- **Day 2:** Multivariate methods for one table (non-supervised / supervised)
- **Day 3:** Multivariate methods to integrate multiple tables (I)
- **Day 4:** Multivariate methods to integrate multiple tables (II)
- **Day 5:** Conditional models, Mendelian Randomization and mediation analysis

## Data

Data for reproducing the slides (folder 'data') and doing the exercises (folder 'data_exercises') can be found [here](https://drive.google.com/file/d/0B9z3DpVIA7h6WEx0Q3Jkd00tYm8/view?usp=sharing). 

Data for final illustration using Markdown can be found [here](https://drive.google.com/file/d/0B9z3DpVIA7h6YW1iOXZnWmZKY0k/view?usp=sharing).

## Example PCA and Hierarchical clustering in diet

Have a look at this paper: [A comparison of the dietary patterns derived by principal component analysis and cluster analysis in older Australians](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4772350/)