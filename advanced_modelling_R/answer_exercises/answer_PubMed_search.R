# Exercise 1
library(RISmed)
copd <- EUtilsSummary("physical activity[ti] + copd[ti] + 
                       Lancet[jo]",
                       db = "pubmed")

copd <- EUtilsSummary("physical activity[ti] + copd[ti]",
                      db = "pubmed")


QueryCount(copd)

# How many last year
metadata.copd <- EUtilsGet(copd)

YearPubmed(metadata.copd)
table(YearPubmed(metadata.copd))

# Get Abstracts

abstracts <- AbstractText(metadata.copd)

abstracts[4]

# How many of these studies are clinical trials
class(abstracts)

ff <- function(x){
  ans <- grep("trial", x)
  lapply(ans, function(x) any(unlist(x)==1))
}

ans <- sapply(abstracts, ff)
sum(unlist(ans))


# Find a paper with a mea-analysis

copd.meta <- EUtilsSummary("physical activity[ti] + copd[ti] + 
                           meta-analysis", db = "pubmed")
QueryCount(copd.meta)
meatadata.copd.meta <- EUtilsGet(copd.meta)

PMID(meatadata.copd.meta)


