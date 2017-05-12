#' Linear model
#'
#' This function estimates model parameters in a linear regression model
#' @param formula formula to be fitted
#' @param data data
#' @param ... other arguments
#' @export

lmModFinal <- function(formula, data, ...){

  mf <- model.frame(formula=formula, data=data)
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, data=mf)
  y <- model.response(mf)

  qx <- qr(x)
  coef <- solve.qr(qx, y, ...)
  df <- nrow(x) - ncol(x)
  sigma2 <- sum((y-x%*%coef)^2)/df
  vcov <- sigma2 * chol2inv(qx$qr)
  colnames(vcov) <- rownames(vcov) <- colnames(x)
  ans <- list(coefficients = coef, vcov = vcov,
              sigma = sqrt(sigma2), df = df,
              data = data.frame(x=x, y=y))
  ans
}
