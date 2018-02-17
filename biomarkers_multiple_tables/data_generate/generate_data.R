# diet.rda copia de base_simulada.dta

library(readstata13)
dd <- read.dta13("c:/juan/cursos/Multi_omic_Madrid/data/base_simulada.dta")
names(dd)

X1 <- dd[,1:16]  # phenotypes and covariates
X2 <- dd[,17:26] # nutrients
X3 <- dd[,27:48] # foods

save(X1, X2, X3, file="c:/Juan/cursos/Multi_omic_Madrid/data_exercises/diet.Rdata")
