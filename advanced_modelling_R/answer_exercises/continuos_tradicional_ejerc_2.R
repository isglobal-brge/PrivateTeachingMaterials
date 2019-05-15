setwd("c:/juan/CREAL/GitHub/PrivateTeachingMaterials/advanced_modelling_R/")

dat <- read.delim("data/agudezavisual.txt")
head(dat)

dat.agg <- aggregate( agudeza ~ indi + luz + group, data=dat, FUN=mean)

dim(dat)
dim(dat.agg)


# Visual inspection

xyplot(agudeza ~ luz|group, group=indi, data=dat.agg,
       type="l", col="blue", lwd=2)


xyplot(agudeza ~ luz|group+ojo, group=indi, data=dat,
       type="l", col="blue", lwd=2)

# ANOVA medidas repetidas con datos agregados

mod1 <- aov(agudeza ~ luz*group + Error(indi), data=dat.agg)
summary(mod1)

# ANOVA medidas repetidas con datos originales

mod2a <- aov(agudeza ~ luz*group + Error(indi), data=dat)
summary(mod2a)

mod2b <- aov(agudeza ~ luz*group*ojo + Error(indi), data=dat)
summary(mod2b)

# Para interpretar miramos los coeficientes 
modLm<- lm(agudeza ~ luz*group*ojo, data=dat)
summary(modLm)

# Plot means

dat.agg2 <- aggregate( agudeza ~ ojo + group, data=dat, FUN=mean)


library(ggplot2)
tgc <- Rmisc::summarySE(dat, measurevar="agudeza", groupvars=c("ojo","group"))
tgc
pd <- position_dodge(0.1) 
ggplot(tgc, aes(x=ojo, y=agudeza, colour=group, group=group)) + 
  geom_errorbar(aes(ymin=agudeza-ci, ymax=agudeza+ci), colour="black", 
                width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=3) 


# ANOVA medidas repetidas con datos originales dos niveles de correlación
#       Error(id/ojo)

mod4 <- aov(agudeza ~ luz*group*ojo + Error(indi/ojo), data=dat)
summary(mod4)
