load("c:/Juan/CREAL/rexposome/nhanes/nh_99-06.Rdata")
table(VarDescription$category)

VarDescription[VarDescription$category=="disease",]
VarDescription[VarDescription$category=="body measures",]
VarDescription[VarDescription$category=="volatile compounds",]

# Data for exercise

nut <- unique(VarDescription[VarDescription$category=="nutrients", "var"])
air <- unique(VarDescription[VarDescription$category=="volatile compounds", "var"])
vars <- c(air, nut, "male", "education",
          "BMXBMI", "RIDAGEYR",
          "colon_cancer_self_report", "any_ht")
bd <- MainTable[,vars]

if(FALSE) {  # set to TRUE to impute data
 set.seed(12345)
 library(mice)
 bd.i <- mice(bd, m=1, method="pmm")
 save(bd.i, file="bd_i.Rdata")
}

load("bd_i.Rdata")

bd.i.c <- complete(bd.i, action=1)
sel <- sapply(bd.i.c,function(x) !any(is.na(x)))
bd.i.c2 <- bd.i.c[, sel]

o <- complete.cases(bd[,87:88])
nhanes <- bd.i.c2[o,]
nn <- names(nhanes)

vars2 <- vars[vars%in%nn]

nhanes.nut <- nhanes[,vars2[!vars2%in%air]]
nhanes.air <- nhanes[,vars2[!vars2%in%nut]]
names(nhanes.nut)[32:35] <- c("BMI", "AGE", "COLON", "HT")
names(nhanes.air)[48:51] <- c("BMI", "AGE", "COLON", "HT")

air.desc <- VarDescription[VarDescription$category=="volatile compounds",
               c("var", "var_desc")]
nut.desc <- VarDescription[VarDescription$category=="nutrients",
                           c("var", "var_desc")]


save(nhanes.nut, nhanes.air, air.desc, nut.desc,
     file="nhanes.Rdata")

load("nhanes.Rdata")
set.seed(12345)
nhanes.nut$CHOL <- with(nhanes.nut, round(170 + 2.3*log(LBXB12) + 5.56*log(LBXFOL)
                        + 2*log(LBXLUZ) -4*log(LBXPHF) - 3.2*log(LBXPHE)
                        + rnorm(nrow(nhanes.nut), 0, 2),0))
summary(nhanes.nut$CHOL)
nhanes.nut$CHOL2 <- cut(nhanes.nut$CHOL, c(-Inf, 200, Inf), labels=c("normal", "hypercol"))

nhanes.air$CHOL <- nhanes.nut$CHOL
nhanes.air$CHOL2 <- nhanes.nut$CHOL2

save(nhanes.nut, nhanes.air, air.desc, nut.desc,
     file="data_exercises/nhanes.Rdata")
