set.seed(123456)

n <- 100

dd <- data.frame(
        age = rnorm(n, 30, 12),
        sex = sample(c("male", "female"), n, rep=TRUE),
        bmi = rnorm(n, 25, 6),
        smoke = sample(c("never", "ex-smoker", "smoker"), n, rep=TRUE),
        cholesterol = rnorm(n, 180, 9)
)

group <- as.factor(sample(c("A", "B"), n, rep=TRUE))


compareGroups <- function(x, y, pvalue=TRUE){
  
  if(missing(x) | missing(y))
    stop("x and y are required")
  
  if(nrow(x)!=length(y))
    stop("x and y should have the same number of observations")
  
  if(!is.factor(y) | length(levels(y))!=2)
    stop("y must be a 2 level factor")
  
  ans <- list()
  
  for (i in 1: ncol(x)){
    print(i)
    xx <- x[,i]
    ans[[i]] <- list()
    ans[[i]]$table <- descriptive(xx, y)
    if (pvalue){
      ans[[i]]$pvalue <- pval(xx, y)
    }
  }
  names(ans) <- colnames(x)
  ans
}


descriptive <- function(x, y) {
  if (is.numeric(x)) {
    m <- aggregate(x ~ y, FUN=mean)
    s <- aggregate(x ~ y, FUN=sd)
    ans <- cbind(m, s[,2])
    colnames(ans) <- c("group", "mean", "sd")
  }
  else{
    tt <- table(y,x)
    pp <- prop.table(tt)
    ans <- cbind(tt, pp)
  }
    ans
}

pval <- function(x, y) {
  if (is.numeric(x)) {
    p <- t.test(x ~ y)$p.value
  }
  else{
    tt <- table(y,x)
    p <- chisq.test(tt)$p.value
  }
  p
}




















































