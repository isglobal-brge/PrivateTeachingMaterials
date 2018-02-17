#
# Task 1
#

library(readstata13)
setwd("c:/juan/cursos/Multi_omic_Madrid/") # change this!!!!
diet <- read.dta13("data_exercises/diet.dta")


# question 1
nrow(diet)

# question 2
diet.m <- subset(diet, estudios=="Bachiller/BUP/COU")
nrow(diet.m)

# queston 3
diet.m[9:12,]

# question 4
median(diet$peso)

# question 5

aggregate(diet$peso, by=list(diet$tipocancer), mean)
aggregate(peso ~ tipocancer, data=diet, mean)

aggregate(diet$peso, by=list(diet$tipocancer, diet$sexo), mean)
aggregate(peso ~ tipocancer + sexo, data=diet, mean)


# question 6
boxplot(t_zinc ~ droplevels(tipocancer), data=diet, xlab="Tipo de cancer", 
        ylab="zinc",  col="red", names=c("Control", "Colorectal", "Breast", 
                                       "Prostate", "Stomach"))
title("Zinc distribution accross tumors")



boxplot(t_zinc ~ droplevels(tipocancer) + sexo, data=diet, xlab="Tipo de cancer", 
        ylab="zinc",  col=c("red", "blue"), las=2)
title("Zinc distribution accross tumors")


#
# TASK 2
#

library(compareGroups)

# question 1
diet.cc <- subset(diet, tipocancer%in%c("Control", "Colorrectal", "Mama"))
table(diet.cc$tipocancer)
dim(diet.cc)

# questin 2
descr1 <- compareGroups(tipocancer ~ edad + sexo + estudios + peso + 
                                    altura + mets_10a + mets_5a + 
                                    Diabetes + Hipertensio + Colesterol, 
                       data = diet.cc )
tab1 <- createTable(descr1)
tab1


# task 3
descr.colon <- compareGroups(tipocancer ~ . -id -casoc -casom -casop -casoe, 
                        data = diet.cc, 
                        subset = tipocancer%in%c("Control", "Colorrectal"))
tab.colon <- createTable(descr.colon, show.ratio=TRUE, 
            show.p.overall=FALSE, show.p.trend = TRUE)
tab.colon
