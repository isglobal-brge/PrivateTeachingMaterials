library(readstata13)
diet <- read.dta13("c:/Juan/CREAL/GitHub/biomarkers_multiple_tables/data/diet2.dta")
names(diet)

diet.cc <- diet[!is.na(diet$casoc),]
dim(diet.cc)
diet.cc$Colesterol <- droplevels(diet.cc$Colesterol)

mod.X <- glm(casoc ~ gra_conveniencefood,data=diet.cc, family="binomial")
summary(mod.X)

mod.M <- glm(Colesterol ~ gra_conveniencefood, data=diet.cc, family="binomial")
summary(mod.M)

mod.Y <- glm(casoc ~ gra_conveniencefood + Colesterol, data=diet.cc, family="binomial")
summary(mod.Y)

mediation::mediate(mod.M, mod.Y, treat = "gra_conveniencefood",
                   mediator = "Colesterol", sims=200)


#
#  Antoher example:  vegetales -> fibra -> cancer
#

diet.cc$t_energy <- scale(diet.cc$t_energy)
diet.cc$gra_veg <- scale(diet.cc$gra_veg)
diet.cc$t_fiber <- scale(diet.cc$t_fiber)

mod.X <- glm(casoc ~ gra_veg + t_energy, data=diet.cc, family="binomial")
summary(mod.X)

mod.M <- glm(t_fiber ~ gra_veg  + t_energy, data=diet.cc, family="gaussian")
summary(mod.M)

mod.Y <- glm(casoc ~ gra_veg + t_fiber + t_energy, data=diet.cc, family="binomial")
summary(mod.Y)

mm <- mediation::mediate(mod.M, mod.Y, treat = "gra_veg",
                   mediator = "t_fiber", sims=500)
summary(mm)
