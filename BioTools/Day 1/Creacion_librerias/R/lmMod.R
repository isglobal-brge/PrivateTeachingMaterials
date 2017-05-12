lmMod <- function(x, ...) UseMethod("lmMod")


# anyadimos un metodo por defecto 'lmMod.default'
lmMod.default <- function(x, y, ...)
{
  # controles !!!

  x <- as.matrix(x)
  y <- as.numeric(y)
  ans <- lmEst(x, y)

  # mas calculos utiles
  ans$fitted.values <- as.vector(x %*% ans$coefficients)
  ans$residuals <- y - ans$fitted.values
  ans$call <- match.call()
  ans$data <- data.frame(y, x)

  # definimos la clase
  class(ans) <- "lmMod"
  ans # devuelve una lista!!!
}


print.lmMod <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
}
