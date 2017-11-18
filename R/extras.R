is.between = function(x,int)  int[1] <= x & x <= int[2]

####
#### The following set of functions is designed to simulate
#### a sample of size n with TAIL distribution 1/(xlog(x))
#### This distribution is used by Embrechts (p.194) to 
#### describe Hill Horror-Plots
####

buildMomentDerConstFunc =  function(D,J1,J2)
{
  output  <- list()
  J <- max(J2)
  J <- if (J == -Inf)  0 else J
  
  if (length(J2)>= 1)
    for (i in 1:length(J2))
      output[[i]] = eval(substitute (function(x,...) x^(J-j)/factorial(D-j),list(j=J2[i], D=D, J=J) ))  
  
  K = length(output)
  
  if (length(J1)>= 1)
    for (i in 1:length(J1))
      output[[K + i]] = eval(substitute (function(x,...) gamma(j+1)/gamma(j+D+1)*x^(j+J),list(j = J1[i], D=D, J=J)))  
  
  
  return(output)
}