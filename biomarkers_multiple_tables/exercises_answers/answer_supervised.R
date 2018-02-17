#
# Task 1
#

# question 1
library(readstata13)
setwd("c:/Juan/CREAL/GitHub/biomarkers_multiple_tables/")
diet <- read.dta13("data_exercises/diet.dta")
vars <- !names(diet)%in%c("id", "casoc", "casom", 
                          "casop", "casoe")
diet <- diet[ , vars]
diet <- diet[complete.cases(diet),]
table(diet$tipocancer)
diet$tipocancer <- droplevels(diet$tipocancer)
table(diet$tipocancer)


# question 2
set.seed(12345)
sel <- sample(1:nrow(diet), 4000)
train <- diet[sel, ]
test <- diet[-sel,]

# question 3
library(randomForest)
mod <- randomForest(tipocancer ~ ., data=train)
predict.rf <- predict(mod, test, type = "class" )
predict.rf
table(obs=test$tipocancer, pred=predict.rf)

flexclust::randIndex(table(obs=test$tipocancer, pred=predict.rf))
