rm(list=ls())

library(magrittr)
load("data/log_Hill_Horror_dist.Rdata")

###
### 2. Set the parameter to vary and create a design of experiment
### Interval over which to maximize
a <- as.numeric(quantile(sample,0.8)) # Threshold a for GLP
m <- 0:3
d <- 3:1

##
## 3. True values of all the parameters at a
##

d1    = DistributionPty::Dlhorror(a, 1)
d2    = -DistributionPty::Dlhorror(a, 2)
d3    = DistributionPty::Dlhorror(a, 3)
m0 = 1-DistributionPty::Dlhorror(a, 0)
m1 = DistributionPty::UpperTrunclHorrorMoment(a,1, n=1e5)
m2 = DistributionPty::UpperTrunclHorrorMoment(a,2, n=1e5)
m3 = DistributionPty::UpperTrunclHorrorMoment(a,3, n=1e5)

truth  = data.frame(d3,-d2,d1,m0,m1,m2,m3)

# Table containing all the possible scenarios when estimating the parameters. 
# For example, the first row of designTab indicates that in the first scenario
# all the derivatives, the tail distribution and the first are estimated.
# In the second scenario, all parameters are estimated except the third derivative.
designTab <- expand.grid(d3 = c(TRUE,FALSE),
                         d2 = c(TRUE,FALSE),
                         d1 = c(TRUE,FALSE),
                         m0 = c(TRUE,FALSE),
                         m1 = c(TRUE,FALSE),
                         m2 = c(TRUE,FALSE),
                         m3 = c(TRUE,FALSE))

# Obtain boostrap sample of each quantity to estimate, i.e.
# F'''(a), F''(a), F'(a), F(a), E[XI(X)>a], E[X^2I(X)=>a],  E[X^3I(X)=>a]
# The moment are obtain via moment estimation, the derivates are obtained via kernel estimation
bootSample <- GLP::getCIMomentAndDerivatives(sample, a,m,d, mc.cores = 4, bootSample = TRUE)$bootSample

# Using the bootstrapped sample, obtain 95% CI and point estimations of the above quantities
# in every possible scenario. A bonferroni correction is applied to account for the number of
# esimated parameters
optim_param <- data.frame(matrix(NA,ncol = 2*length(c(m,d))+1,nrow = 2^length(c(m,d))  ))
names(optim_param) <- c(paste(rep(names(designTab[,1:length(c(d,m))]),each = 2),c("L","U"),sep=""),"cover")

# The level of confidence for the joint interval is 1-alpha
alpha <- 0.05

# Compute the confidence intervals
for (i in 1:nrow(designTab)){
  
  index <- which(designTab[i,] %>% as.matrix)
  
  if (length(index) != 0)
  {
    # Bonferroni correction factor
    bonf <- length(index)
    
    param <- apply(bootSample, 2, quantile, probs = c(alpha/(2*bonf),1-alpha/(2*bonf))) 
    optim_param[i, c(2*index-1, 2*index)] <- as.numeric(c(param[1, index], param[2, index]))
    
    # Does the 1-alpha CI covers the true value?
    optim_param[i, "cover"] = all((param[1,] <= truth & truth <=  param[2,])[,index])
    
    # If the tail distribution is not estimation, set its bounds to 0 and 1
    if (!(4 %in% index))  optim_param[i, c("m0L","m0U")] <- c(0,1)
  }else{
    optim_param[i, c("m0L","m0U")] <- c(0,1)
    optim_param[i, "cover"] = TRUE
  }
}

save(bootSample,optim_param,truth, designTab, file = "data/log_Hill_data_optimparam.Rdata")
