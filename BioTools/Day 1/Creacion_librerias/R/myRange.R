myrange<-function(x)
{
  x<-sort(x)
  ans<-x[length(x)]-x[1]  
  out<-list(x=x, rango=ans)  
  class(out)<-"myrange"
  out
}

print.myrange<-function(x)
{
  cat("El rango de tus valores es: \n")
  cat(x$rango)
  cat("\n")
}