setwd("c:/juan/CREAL/GitHub/PrivateTeachingMaterials/advanced_modelling_R/")

dat <- read.delim("data/agudezavisual.txt")
head(dat)

library(nlme)
dat.s <- groupedData(agudeza ~ luz | ojo/indi,
                     data=dat)
dat.s <- groupedData(agudeza ~ luz | indi,
                     data=dat)


plot(dat.s)
plot(dat.s, inner="group")

library(lattice)
xyplot(agudeza ~ luz | group, groups = indi, data=dat,
       type="l", col="darkgray", lwd=2)

mod.int <- lme(agudeza ~ luz*group, dat.s, random = ~1)
summary(mod.int)

mod.slope <- lme(agudeza ~ luz*group, dat.s, random = ~ luz)
summary(mod.slope)


plot(augPred(mod.int, level = 0:1, length.out = 2))



### ejercicio 2

data("Milk", package="nlme")
head(Milk)
dim(Milk)


df <- as.data.frame(Milk)
df.agg <- aggregate(protein ~ Time + Diet, df, mean)

ggplot(df, aes(x = Time, y = protein, color=Diet)) +
  geom_line(aes(group = Cow), alpha = .3) +
  geom_line(data = dd.agg, alpha = .8, size = 3) +
  scale_color_manual(values=c("darkgreen", "green", "lightgreen"))

Milk$time2 <- cut(Milk$Time, c(0, 4, Inf), 
                  labels=c("Lactancia", "No lactancia"))

mod <- lme(protein ~ time2*Diet, data = Milk, random = ~ 1)
summary(mod)


mod2 <- lme(protein ~ Time*Diet, data = Milk, random = ~ 1)
summary(mod2)

BIC(mod)
BIC(mod2)

