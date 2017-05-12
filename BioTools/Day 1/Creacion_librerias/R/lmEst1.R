lmEst1 <- function(x, y, method = "OLS"){
  # control errores
  method.type <- c("OLS", "lm")
  m <- charmatch(method, method.type, nomatch = 0)
  if (m==0)
    stop("method should be 'OLS' or 'lm' ")
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