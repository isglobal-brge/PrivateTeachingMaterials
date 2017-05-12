ec_grado2 <- function(a, b, c)
{
  # Calculamos el discriminante (b^2 - 4ac)  disc <- (b^2)-(4*a*c)
  # Primera solución
  sol1 <- (-b + sqrt(disc))/(2*a)
  # Segunda solución
  sol2 <- (-b - sqrt(disc))/(2*a)
  # Devolvemos un vector con las dos soluciones y cerramos la llave para indicar el final del grupo de
  # instrucciones que forman la función  
  ans <- c(sol1,sol2)
  ans
}
