lmEst4 <- function(formula, data, method = "OLS"){
  # control errores
  method.type <- c("OLS", "lm")
  m <- charmatch(method, method.type, nomatch = 0)
  if (m==0) {
    warning(paste("method=", method, "is not.supported. Using 'OLS' ..."))
    m <- 1
  }
  
  mf <- model.frame(formula=formula, data=data)
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, data=mf)
  y <- model.response(mf)

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