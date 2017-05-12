lmEst5 <- function(formula, data, contrast = NULL, ...) {

  mf <- model.frame(formula=formula, data=data)
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, data=mf)
  y <- model.response(mf)
  
  ## calcula la descomposicion QR de x
  qx <- qr(x)

  ## calcula (x'x)^(-1)x'y
  n <- as.integer(nrow(qx$qr)) 
  k <- as.integer(qx$rank)
  ny <-as.integer(1)
  
  z <- .Fortran("dqrcf", as.double(qx$qr), n, k, 
                as.double(qx$qraux),  
                as.double(matrix(y, nrow=n, ncol=ny)),
                as.integer(ny), 
                coef =matrix(0,  nrow = k, ncol = ny),
                info = as.integer(0),
                NAOK = TRUE, PACKAGE = "base")
  
  coef <- z$coef
  info <- z$info
  
  df <- nrow(x) - ncol(x)
  sigma2 <- sum((y-x%*%coef)^2)/df
  vcov <- sigma2 * chol2inv(qx$qr)
  colnames(vcov) <- rownames(vcov) <- colnames(x)

  
  # resultado
  ans <- list(coefficients = coef, vcov = vcov, 
              sigma = sqrt(sigma2), df = df, info=info)
  ans
}
