myTest<-function(x,...)  UseMethod("myTest")

myTest.default <- function(x, y, ...){
  #	??? algunos controles  
  xx<-table(x,y)  
  myTest.table(xx, ...)
}

myTest.table<-function(x,...) {
  #	??? algunos controles 
  myTest.matrix(x, ...)
}

myTest.matrix<-function(x, ...) {
  # ??? programa con todos los c??lculos
  # implementa la teor??a ???  
  ans<-fisher.test(x, ...)  
  class(ans)<-"myTest"  
  ans
}

print.myTest<-function(x,...)
{
  cat("Results for my new test \n")  
  cat(" p-value of 'new test':", x$p.value)
  cat("\n")
}
