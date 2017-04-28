is.between = function(x,int)  int[1] <= x & x <= int[2]

####
#### The following set of functions is designed to simulate
#### a sample of size n with TAIL distribution 1/(xlog(x))
#### This distribution is used by Embrechts (p.194) to 
#### describe Hill Horror-Plots
####

buildMomentDerConstFunc =  function(D,d,m)
{
  output  =  list()
  
  if (length(d)>= 1)
    for (j in 1:length(d))
      output[[j]] = eval(substitute (function(x,...) (-1)^(sgn+1)*x^(order)/factorial(order), list(sgn = d[j], order=D - d[j]) ))  
    
    K = length(output)
    
    if (!(0 %in% m)){ 
      output[[K + 1]] = eval(substitute (function(x,...) x^order/factorial(order), list(order=D) ))
      K = K+1
    }

      if (length(m)>= 1)
        for (j in 1:length(m))
          output[[K + j]] = eval(substitute (function(x,...) factorial(order)/factorial(D+order)*x^(order+D),
                                          list(order = m[j])))  
    
    
    return(output)
}