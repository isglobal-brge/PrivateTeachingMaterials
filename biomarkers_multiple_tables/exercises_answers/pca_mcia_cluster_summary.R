library(made4)
library(omicade4)
library(readstata13)
dd <- read.dta13("c:/Juan/CREAL/GitHub/biomarkers_multiple_tables/data/diet2.dta")
names(dd)

# PCA 
dd.ok <- t(scale(dd[,18:50]))
mm <- ord(dd.ok)
plotgenes(mm)
topgenes(mm, axis=1, end="pos")

# CIA
dd.ok2 <- list(nut=dd.ok[1:10,], food=dd.ok[11:33,])
mm2 <- mcia(dd.ok2)
plot(mm2, df.color = c("red", "blue"), gene.nlab = 5)
plotVar(mm2, nlab = 5)

# Hierarchical clustering
plot(hclust(dist(dd.ok)))

