# Prepare system

rtools <- "C:\\Rtools\\bin"
gcc <- "C:\\Rtools\\min_gw_64\\bin" 
path <- strsplit(Sys.getenv("PATH"), ";")[[1]] 
new_path <- c(rtools, gcc, path) 
new_path <- new_path[!duplicated(tolower(new_path))] 
Sys.setenv(PATH = paste(new_path, collapse = ";"))

setwd("src")
system("R CMD SHLIB -o prueba.dll prueba.f")



# load data
bass <- read.delim("data/bass.txt")

# data analyses
source("R/lmEst.R")
lmEst(cbind(1, bass$Alkalinity), bass$Mercury)
lm(Mercury ~ Alkalinity, data=bass)

source("R/lmEst1.R")
lmEst1(cbind(1, bass$Alkalinity), bass$Mercury, method="qr")

source("R/lmEst2.R")
lmEst2(cbind(1, bass$Alkalinity), bass$Mercury, method="qr")

source("R/lmEst3.R")
lmEst3(cbind(1, bass$Alkalinity))
lmEst3(cbind(1, bass$Alkalinity), bass$Mercury[-1])
lmEst3(cbind(1, bass$Alkalinity), bass$Mercury)

lmEst(cbind(1, bass$Alkalinity), bass$Mercury, tol=1e-8)

source("R/lmEst4.R")
lmEst4(Mercury ~ Alkalinity, bass)
lmEst4(Mercury ~ Alkalinity - 1, bass)
lmEst4(Mercury ~ pH, bass)
lmEst4(Mercury ~ . , bass[,-1])

source("R/lmEst5.R")
is.loaded("dqrcf")  
lmEst5(Mercury ~ Alkalinity, bass)

# dll

system("R CMD SHLIB -o prog.dll prog.f") 


# classes and methods

set.seed(123456)
casco <- sample(c("ca","co"), 200, replace=TRUE)
fuma <- sample(c("si","no"), 200, replace=TRUE)
tt<-table(casco, fuma)
tt

source("R/myRange.R")
myrange(c(1,5,6,9,2,19,3))


source("R/myTest.R")
myTest(casco, fuma)
myTest(tt)
mm<-matrix(tt, nrow=2, ncol=2)
myTest(mm)

# let's turn on lineal model (after defining the .default)
source("R/lmMod.R")
source("R/plotlmMod.R")
lmMod(cbind(1, bass$Alkalinity), bass$Mercury)

mod<-lmMod(cbind(Const=1, Alk=bass$Alkalinity),  bass$Mercury)
coef(mod)
fitted(mod)
resid(mod)
plot(mod)


# Creating libraries and manual with devtools and roxygen2
library(devtools)
library(roxygen2)
setwd("C:/Juan/CREAL/cursos/BioTools/Day_1")
create("lmPackage")

# Pongo mis funciones en la carpeta 'R'
# Pongo los datos en 'data' o 'extdata' si son ficheros de texto

setwd("lmPackage")
document()





