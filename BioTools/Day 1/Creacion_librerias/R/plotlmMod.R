plot.lmMod <- function(x, ...)
{
  data <- x$data  
  coefs <- coef(x)
  lab <- colnames(data)
  plot(data[,3], data[,1], xlab=lab[3], ylab=lab[1], ...)  
  abline(coefs[1], coefs[2], lwd=2)
  ff <- paste(lab[1], "=", round(coefs[1],2) , "+",
              round(coefs[2],2), lab[3], sep="")
  legend("topright", ff, bty="n", cex=0.8)
}