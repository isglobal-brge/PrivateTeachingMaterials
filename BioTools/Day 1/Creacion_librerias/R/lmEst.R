lmEst <- function(x, y) {
  ## calcula la descomposici??n QR de x
  qx <- qr(x)

  ## calcula (x'x)^(-1)x'y
  coef <- solve.qr(qx, y)
  
  #grados de libertad y desviaci??n standard de los residuales
  df <- nrow(x) - ncol(x)
  sigma2 <- sum((y-x%*%coef)^2)/df
  
  ## calcula sigma^2 * (x???x)^-1
  vcov <- sigma2 * chol2inv(qx$qr)
  colnames(vcov) <- rownames(vcov) <- colnames(x)
  
  # resultado
  ans <- list(coefficients = coef, vcov = vcov, 
              sigma = sqrt(sigma2), df = df)
  ans
}
