rm(list=ls())

RutinesLocals<-"C:/programs/Dropbox/rutines"

source(file.path(RutinesLocals,"carrega.llibreria.r"))
source(file.path(RutinesLocals,"table2.r"))
source(file.path(RutinesLocals,"subset2.r"))
source(file.path(RutinesLocals,"merge2.r"))
source(file.path(RutinesLocals,"add.cases.r"))
source(file.path(RutinesLocals,"format2.r"))
source(file.path(RutinesLocals,"read.spss4.r"))
source(file.path(RutinesLocals,"spss_varlist.r"))
source(file.path(RutinesLocals,"intervals.r"))
source(file.path(RutinesLocals,"prepare.r"))
source(file.path(RutinesLocals,"export.SPSS.r"))

library(compareGroups)

library(xlsx)

library(gdata)


setwd("C:/programs/Dropbox/cursShiny/data")

data(predimed)
data(regicor)
data(SNPs)


################ PREDIMED #############

set.seed(123456)
datos <- predimed[sample(1:nrow(predimed), 200),]


### TXT ###

write.table(datos, file="predimed.dat", sep="\t", dec=".", row.names=FALSE)

write.csv(datos, file="predimed.csv", row.names=FALSE)



### EXCEL (2003) ###

write.xlsx(datos, "predimed.xlsx", row.names = FALSE, showNA = FALSE)

# "predimed.xls"  fer-ho des d'Excel 


### SPSS ###

temp<-predimed
for (j in 1:ncol(temp)){
  attr(temp[,j],"vari.label")<-Hmisc::label(temp[,j])
  if (is.factor(temp[,j])){
    ll<-levels(temp[,j])
    temp[,j]<-as.integer(temp[,j])
    attr(temp[,j],"value.labels")<-structure(1:length(ll),names=ll)
  }
}

export.SPSS(temp, file.save = file.path(getwd(),"predimed.sav"), file.runsyntax = "C:/programs/PSPP/bin/pspp.exe") 



################ REGICOR #############

set.seed(123456)
datos <- regicor[sample(1:nrow(regicor), 200),]


### TXT ###

write.table(datos, file="regicor.dat", sep="\t", dec=".", row.names=FALSE)

write.csv(datos, file="regicor.csv", row.names=FALSE)



### EXCEL (2003) ###

write.xlsx(datos, "regicor.xlsx", row.names = FALSE, showNA = FALSE)

# "regicor.xls"  fer-ho des d'Excel 

### SPSS ###

temp<-regicor
for (j in 1:ncol(temp)){
  attr(temp[,j],"vari.label")<-label(temp[,j])
  if (is.factor(temp[,j])){
    ll<-levels(temp[,j])
    temp[,j]<-as.integer(temp[,j])
    attr(temp[,j],"value.labels")<-structure(1:length(ll),names=ll)
  }
}

export.SPSS(temp, file.save = file.path(getwd(),"regicor.sav"), file.runsyntax = "C:/programs/PSPP/bin/pspp.exe") 



#### multiple sheets
data(regicor)
data(predimed)
data(SNPs)

set.seed(123456)
regicor <- regicor[sample(1:nrow(regicor), 200),]
predimed <- predimed[sample(1:nrow(predimed), 200),]
SNPs <- SNPs

write.xlsx(regicor, "datos.xlsx", sheetName="regicor", row.names = FALSE, showNA = FALSE)
write.xlsx(predimed, "datos.xlsx", sheetName="predimed", row.names = FALSE, showNA = FALSE, append = TRUE)
write.xlsx(SNPs, "datos.xlsx", sheetName="SNPs", row.names = FALSE, showNA = FALSE, append = TRUE)




