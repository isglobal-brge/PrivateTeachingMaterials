## ----data2---------------------------------------------------------------
library(pdfCluster)
data(oliveoil)
head(oliveoil)

## ----dataTT--------------------------------------------------------------
set.seed(1234)
ss <- sample(1:nrow(oliveoil), 200)
train <- oliveoil[-ss,-2]
test <- oliveoil[ss,-2]

## ----lda-----------------------------------------------------------------
library(MASS)
olive.lda <- lda(macro.area~., train)
pregion.lda <- predict(olive.lda, test)$class
table(test[,1], pregion.lda)

## ----ldaPlot-------------------------------------------------------------
plot(predict(olive.lda, test)$x,
     col=as.numeric(test[,1]))

## ----trees---------------------------------------------------------------
library(rpart)
olive.rp <- rpart(macro.area~., train, 
          method="class")
olive.rp

## ----treesPlot-----------------------------------------------------------
plot(olive.rp)
text(olive.rp)

## ----ldaPredict----------------------------------------------------------
temp <- predict(olive.rp, test)
head(temp)
pregion.rp <- apply(temp, 1, function(x) which(x==1))

## ----ldaPredict2---------------------------------------------------------
table(test[,1], pregion.rp)

## ----svm-----------------------------------------------------------------
library(e1071)
olive.svm <-  svm(macro.area ~. , data = train)
pregion.svm <- predict(olive.svm, test)
table(test[,1], pregion.svm)

## ----nnet----------------------------------------------------------------
library(nnet)
olive.nnet <-  nnet(macro.area ~. , data = train,
               size=2)
pregion.nnet <- predict(olive.nnet, test, type="class")
table(test[,1], pregion.nnet)

## ----nnet2---------------------------------------------------------------
olive.nnet <-  nnet(macro.area ~. , data = train,
               size=4)
pregion.nnet <- predict(olive.nnet, test, type="class")
table(test[,1], pregion.nnet)

## ----boost---------------------------------------------------------------
library(adabag)
olive.boost <-  boosting(macro.area ~. , data = train,
                control = rpart.control(maxdepth = 2))
pregion.boost <- predict(olive.boost, test, type="class")$class
table(test[,1], pregion.boost)

## ----rf------------------------------------------------------------------
library(randomForest)
olive.rf <-  randomForest(macro.area ~. , data = train)
pregion.rf <- predict(olive.rf, test, type="class")
table(test[,1], pregion.rf)

## ----perf----------------------------------------------------------------
library(flexclust)
randIndex(table(test[,1], pregion.rf))
randIndex(table(test[,1], pregion.lda))
randIndex(table(test[,1], pregion.rp))
randIndex(table(test[,1], pregion.boost))

## ----ROC_1---------------------------------------------------------------
library(pROC)
data(aSAH)
head(aSAH)

## ----ROC_2---------------------------------------------------------------
rocobj <- plot.roc(aSAH$outcome, aSAH$s100b,
                main="Confidence intervals", 
                percent=TRUE,
                ci=TRUE,
                print.auc=TRUE)

## ----ROC_3---------------------------------------------------------------
rocobj <- plot.roc(aSAH$outcome, aSAH$s100b,
                main="Confidence intervals", percent=TRUE,
                ci=TRUE, 
                print.auc=TRUE)
ciobj <- ci.se(rocobj, progress = "none",
               specificities=seq(0, 100, 5)) #This can be selected (grid for computing bands)
plot(ciobj, type="shape", col="lightblue") # plot as a blue shape

## ----compareAUC----------------------------------------------------------
rocobj1 <- plot.roc(aSAH$outcome, aSAH$s100,
                    main="Statistical comparison", percent=TRUE, col="blue")
rocobj2 <- lines.roc(aSAH$outcome, aSAH$ndka, percent=TRUE, col="red")
testobj <- roc.test(rocobj1, rocobj2)
text(50, 50, labels=paste("p-value =", format.pval(testobj$p.value)), adj=c(0, .5))
legend("bottomright", legend=c("S100B", "NDKA"), col=c("blue", "red"), lwd=2)

