X <- matrix(rnorm(1000000), nrow=100, ncol=10000)
y <- sample(c("A", "B"), 100, rep=TRUE)

ff <- function(i, datos, y){
  xx <- datos[,i]
  cov <- sd(xx)
  mm <- mean(xx)
  ss <- lm(xx ~ y)
  ans <- list(cov=cov, mean=mm, lm=ss)
  ans
}

out <- mclapply(1:ncol(X), ff, datos=X, y=y)

out <- list()
for (i in 1: ncol(x)){
  xx <- X[,i]
  cov <- sd(xx)
  mm <- mean(xx)
  ss <- lm(xx ~ y)
  ans <- list(cov=cov, mean=mm, lm=ss)
  out[[i]] <- ans
  
}
