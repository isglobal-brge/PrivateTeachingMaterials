## ------------------------------------------------------------------------
library(foreign)
library(Hmisc)
library(readxl)

## ------------------------------------------------------------------------
df<-read.table("data/parto2.dat", sep=";", as.is=TRUE, header=FALSE)
head(df)

## ------------------------------------------------------------------------
df<-read_excel("data/mujeres.xlsx")
class(df)
class(df) <- "data.frame"
head(df)

## ------------------------------------------------------------------------
df <- read.dta("data/partoFin.dta",
             convert.dates = TRUE, convert.factors = TRUE)
head(df)

## ---- eval=FALSE---------------------------------------------------------
## library(readstata13)
## df <- read.dta13("data/partoFin.dta",
##              convert.dates = TRUE, convert.factors = TRUE)

## ------------------------------------------------------------------------
df <- spss.get("data/parto2.sav",use.value.labels=TRUE, allow="_", 
             datevars=c("dia_nac","dia_entr","ulti_lac"))
head(df)

## ---- eval=FALSE---------------------------------------------------------
## write.table(df,"parto2ex.dat")

## ---- eval=FALSE---------------------------------------------------------
## write.dta(df, file="c:/juan/data/bd.dta"), version=7L)
## save.dta13(df, file="c:/juan/data/bd.dta")

## ---- eval=FALSE---------------------------------------------------------
## save(df, file="c:/juan/data/bd.Rdata")) # or .rda

## ---- eval=FALSE---------------------------------------------------------
## load("c:/juan/data/bd.Rdata")) # an object df will be in R

## ------------------------------------------------------------------------
library(Hmisc)
df <- spss.get("data/partoFin.sav", allow="_", 
               datevars=c("dia_nac", "dia_entr", "ulti_lac"))

## ------------------------------------------------------------------------
head(df)

## ------------------------------------------------------------------------
nrow(df)
ncol(df)

## ------------------------------------------------------------------------
names(df)

## ------------------------------------------------------------------------
summary(df)

## ------------------------------------------------------------------------
df$sexo

## ------------------------------------------------------------------------
df[,2]

## ------------------------------------------------------------------------
df[,c("sexo", "peso", "edad")]

## ------------------------------------------------------------------------
df[,c(1,3,5)]

## ------------------------------------------------------------------------
df[4,]

## ------------------------------------------------------------------------
df[4:10,]

## ------------------------------------------------------------------------
subset(df, sexo=="niña")

## ------------------------------------------------------------------------
table(df$naci_ca)
subset(df, naci_ca%in%c("Española", "Otras"))

## ------------------------------------------------------------------------
mean(df$edad)

## ------------------------------------------------------------------------
sd(df$edad)

## ------------------------------------------------------------------------
median(df$edad)

## ------------------------------------------------------------------------
quantile(df$edad, c(0.25, 0.50, 0.75))

## ------------------------------------------------------------------------
with(df, cor(peso, edad))

## ------------------------------------------------------------------------
with(df, cor(peso, edad, method="spearman"))

## ------------------------------------------------------------------------
hist(df$peso)

## ------------------------------------------------------------------------
plot(df$sexo)

## ------------------------------------------------------------------------
boxplot(df$peso, ylab="Peso (kgs.)")

## ------------------------------------------------------------------------
boxplot(peso ~ sexo , data=df, col="red", 
        ylab="Peso (kgs.)", xlab="Sexo")

## ------------------------------------------------------------------------
plot(peso ~ edad, data=df, col=sexo, pch=19)
title("Weight by mother age")
legend("topright", c("boy","girl"), fill=c(1,2))

## ------------------------------------------------------------------------
t.test(df$peso, mu=4)

## ------------------------------------------------------------------------
t.test(peso ~ sexo, data=df)

## ------------------------------------------------------------------------
t.test(df$horas_an, df$horas_de, paired = TRUE)

## ------------------------------------------------------------------------
freq <- with(df, table(sexo, tip_par))
fisher.test(freq)

## ------------------------------------------------------------------------
cor.test(df$peso, df$edad)

## ------------------------------------------------------------------------
cor.test(df$peso, df$edad, method="spearman")

## ------------------------------------------------------------------------
model <- lm(peso ~ edad, data=df)
summary(model)

## ------------------------------------------------------------------------
model <- glm(tip_par ~ edad, data=df, family="binomial")
summary(model)

## ---- scripting----------------------------------------------------------
sel <- df$tx=="Intensivo"
sel[1:6]
df.intensive <- df[sel,]
model.int <- glm(tip_par ~ edad, data=df.intensive, 
                 family="binomial")
summary(model.int)

## ------------------------------------------------------------------------
load("data/russett.Rdata")
head(X_agric)
head(X_ind)

## ------------------------------------------------------------------------
head(X_polit)
X <- list(tab1 = X_agric, tab2 = X_ind, tab3 = X_polit)
length(X)
head(X[[1]])

## ---- results="hide"-----------------------------------------------------
library(compareGroups)
data(predimed)
# ?predimed

## ------------------------------------------------------------------------
head(predimed)

## ------------------------------------------------------------------------
descr <- compareGroups(group ~ sex + age + smoke, predimed)
descr

## ------------------------------------------------------------------------
descr <- compareGroups(group ~ ., predimed)
descr

## ------------------------------------------------------------------------
descr2 <- compareGroups(group ~ . -sex -age -event, predimed)
descr2

## ------------------------------------------------------------------------
descrtable <- createTable(descr)
descrtable

## ------------------------------------------------------------------------
update(descrtable, hide.no='no')

## ------------------------------------------------------------------------
update(descrtable, hide.no='no', show.n = TRUE)

## ------------------------------------------------------------------------
update(descrtable, hide.no='no', show.n = TRUE, type=1)

## ------------------------------------------------------------------------
descr <- update(descr, method=2)
createTable(descr, hide.no="no")

## ------------------------------------------------------------------------
update(descrtable, hide.no='no', digits=1, digits.p=5)

## ------------------------------------------------------------------------
descr <- update(descr, method=c(age=2, p14=2))
createTable(descr, hide.no="no")

## ------------------------------------------------------------------------
table(predimed$hyperchol)

## ------------------------------------------------------------------------
descr <- compareGroups(hyperchol ~ ., predimed)
createTable(descr, hide.no="no", show.ratio=TRUE, 
            show.p.overall=FALSE, show.p.trend = FALSE)

## ------------------------------------------------------------------------
predimed$tevent <- with(predimed, Surv(toevent, event=="Yes"))

## ------------------------------------------------------------------------
descr <- compareGroups(tevent ~ .- toevent-event, predimed)
createTable(descr, hide.no="no", show.ratio=TRUE, 
            show.p.overall=FALSE)

## ------------------------------------------------------------------------
label(predimed$age) <- "Age of participant"

## ------------------------------------------------------------------------
descrtable <- createTable(compareGroups(group ~ ., predimed))
varinfo(descrtable)

## ------------------------------------------------------------------------
descrtable <- createTable(compareGroups(group ~ ., predimed))
descrtable[c('age','bmi')]
descrtable[c(1,4)]

## ------------------------------------------------------------------------
descr <- compareGroups(group ~ ., predimed)
plot(descr['bmi'])

## ------------------------------------------------------------------------
plot(descr['bmi'], bivar=TRUE)

## ------------------------------------------------------------------------
plot(descr['sex'])

## ------------------------------------------------------------------------
plot(descr['sex'], bivar=TRUE)

## ---- eval=FALSE---------------------------------------------------------
## # CSV
## export2csv(descrtable, file="tabla.csv", sep=";")
## # Excel
## export2xls(descrtable, file="tabla.xlsx")
## # Word
## export2word(descrtable, file="tabla.docx")
## # Latex
## export2tex(descrtable, file="tabla.tex")

## ------------------------------------------------------------------------
export2md(descrtable)

## ------------------------------------------------------------------------
sessionInfo()

