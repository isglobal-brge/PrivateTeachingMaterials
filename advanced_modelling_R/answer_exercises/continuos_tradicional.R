
setwd("c:/juan/CREAL/GitHub/PrivateTeachingMaterials/advanced_modelling_R/")

dat <- read.table("data/recuperainfarto.txt")
head(dat)
names(dat) <- c("grupo", paste0("time", 1:5))
head(dat)
dat$id <- 1:nrow(dat)
head(dat)


# crear datos largos

dat.long <- reshape(dat, direction = "long",
                    varying = list(paste0("time", 1:5)),
                    idvar = c("id","grupo"),
                    v.names="bartel",
                    timevar = "time", times=1:5)
head(dat.long)

# grafico
library(lattice)
xyplot(bartel ~ time | grupo, group=id,
       data = dat.long, type="l", 
       col=ifelse(dat.long$grupo=="B", "gray70", "gray20"),
       lwd=2)

# End-point 
mod1 <- aov(time5 ~  time1 + grupo, data=dat)
summary(mod1)

# MANOVA
mod2 <- manova(cbind(time1, time2, time3, time4, time5) ~ grupo, 
               data=dat)
summary(mod2)


datM <- dat[ , grep("^time", names(dat))]
mod3 <- manova(datM ~ grupo, data=dat)
summary(mod3)


# ANOVA con medidas repetidas

mod4 <- aov(bartel ~ time * grupo + Error(id), data=dat.long)
mod4
summary(mod4)


