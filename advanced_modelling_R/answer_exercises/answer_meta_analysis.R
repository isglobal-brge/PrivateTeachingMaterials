library(meta)

library(readxl)
example <- read_excel("C:/Juan/CREAL/GitHub/PrivateTeachingMaterials/advanced_modelling_R/data/example.xlsx")
example

library(meta)
res <- metacont(Ne, Me, Se, 
                Nc, Mc, Sc, 
                data=example, studlab = Study)

res

summary(res)

### NEXT EXERCISES

res.rmle <- metacont(Ne, Me, Se, 
                    Nc, Mc, Sc, 
                    data=example, studlab = Study,
                    method.tau="REML")

# This makes the same!
res.rmle <- update(res, method.tau="REML")
res.rmle

res.HE <- update(res, method.tau="HE")
res.HS <- update(res, method.tau="HS")

res$tau
res.rmle$tau
res.HE$tau
res.HS$tau

forest(res.rmle, layout = "JAMA")
forest(res.HS, layout = "JAMA")

forest(res.HS, prediction = TRUE)


# NEXT EXERCISE
or <- c(0.12, 0.11, 0.52, 0.49, 0.74, 0.35, 0.19, 0.25)
ciInf <- c(0.01, 0.04, 0.4, 0.29, 0.22, 0.16, 0.06, 0.11)

beta <- log(or)
se <- (beta - log(ciInf))/1.96

res.ex <- metagen(beta, se, sm="OR")
forest(res.ex, prediction = TRUE)

forest(res.ex, col.random = "tomato", col.fixed = "red",
       col.study = c(rep("blue", 2), "red", rep("blue", 5)),
       col.square = c(rep("blue", 2), "red", rep("blue", 5)),
       layout="JAMA")


# NEXT EXERCISE (with 'example' dataset  )

library(altmeta)

res <- metacont(Ne, Me, Se, 
                Nc, Mc, Sc, 
                data=example, studlab = Study)

forest(res)
metainf(res)


outlier <- metaoutliers(res$TE, res$seTE^2)
outlier$outliers

new.res <- metacont(Ne, Me, Se, 
                Nc, Mc, Sc, 
                data=example[-3, ], studlab = Study)


# better 
new.res <- metacont(Ne, Me, Se, 
                    Nc, Mc, Sc, 
                    data=example[-outlier$outliers, ],
                    studlab = Study)

# even better
new.res <- metacont(Ne, Me, Se, 
                    Nc, Mc, Sc, 
                    data = example,
                    subset = Study!=Study[outlier$outliers],
                    studlab = Study)


new.res


# Exclude run meta-analysis exluding this study but keeping the output

new2.res <- metacont(Ne, Me, Se, 
                    Nc, Mc, Sc, 
                    data = example,
                    exclude = Study==Study[outlier$outliers],
                    studlab = Study)


new2.res
forest(new2.res)


# Fleiss example

data("Fleiss93", package = "meta")
Fleiss93

res.fleiss <- metabin(event.e, n.e, event.c, n.c,
                      data=Fleiss93, studlab = study,
                      sm="OR")

res.fleiss

forest(res.fleiss, layout = "JAMA")

metainf(res.fleiss)

out <- metaoutliers(res.fleiss$TE, res.fleiss$seTE^2)
out$outliers
