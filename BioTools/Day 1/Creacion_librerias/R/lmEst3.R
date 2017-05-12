lmEst3 <- function(x, y, method = "OLS"){
  # control errores
  method.type <- c("OLS", "lm")
  m <- charmatch(method, method.type, nomatch = 0)
  if (m==0) {
    warning(paste("method=", method, "is not.supported. Using 'OLS' ..."))
    m <- 1
  }
  
  if (missing(x) | missing(y))
    stop("x and y arguments are required")
  if(nrow(x)!=length(y))
    stop("x and y should have the same length")
  if (any(is.na(x)) | any(is.na(y))) {
    warning("There are missing values. Only complete cases have been analyzed")
    o <- complete.cases(x,y)
    x <- x[o,]
    y <- y[o,]
  }
  
  if (m ==1) {
    # calculos de antes
    qx <- qr(x)
    coef <- solve.qr(qx, y)
    df <- nrow(x) - ncol(x)
    sigma2 <- sum((y-x%*%coef)^2)/df
    vcov <- sigma2 * chol2inv(qx$qr)
    colnames(vcov) <- rownames(vcov) <- colnames(x)
    ans <- list(coefficients = coef, vcov = vcov, 
               sigma = sqrt(sigma2), df = df)
  }
  else { # if (m==2) 
    ans <- lm(x ~ y)
  } 
  ans
}