rm(list=ls())

library(DistributionPty)
library(GLP)
library(magrittr)
setwd("/Users/cmottet/Projects/R/MonotoneDerivatives")
load("data/Hill_Horror_dist.Rdata")
Tsample <- log(sample)

###
### 2. Set the parameter to vary and create a design of experiment
### Interval over which to maximize
a <- as.numeric(quantile(Tsample,0.8)) # Threshold a for GLP
P <- 1- c(10^-1,10^-2,10^-3,10^-4)
q <- qlhorror(P)
m <- 0:3
d <- 3:1
D <- 0:5
direction = 1#<- c(-1,1)

##
## 3. True values of all the parameters at a
##

d1    = Dlhorror(a, 1)
d2    = -Dlhorror(a, 2)
d3    = Dlhorror(a, 3)
m0 = 1-Dlhorror(a, 0)
m1 = UpperTrunclHorrorMoment(a,1, n=1e5)
m2 = UpperTrunclHorrorMoment(a,2, n=1e5)
m3 = UpperTrunclHorrorMoment(a,3, n=1e5)

truth  = data.frame(d3,-d2,d1,m0,m1,m2,m3)

designTab <- expand.grid(d3 = c(TRUE,FALSE),
                        d2 = c(TRUE,FALSE),
                        d1 = c(TRUE,FALSE),
                        m0 = c(TRUE,FALSE),
                        m1 = c(TRUE,FALSE),
                        m2 = c(TRUE,FALSE),
                        m3 = c(TRUE,FALSE),
                        D  = D,
                        "direction" = 1,
                        P  = P)
# As D needs to be larger than the largest derivative,
# we drop the rows that do not meet this constraint
dMax <- apply(designTab[,1:length(d)],1,function(x)max(d[x]))
dMax[dMax == -Inf] <- 0
designTab <- subset(designTab, dMax <= designTab$D) ; rownames(designTab) <- 1:nrow(designTab)


bootSample <- getCIMomentAndDerivatives(Tsample, a,m,d, mc.cores = 2, bootSample = TRUE)$bootSample


optim_param <- data.frame(matrix(NA,ncol = 2*length(c(m,d))+1,nrow = 2^length(c(m,d))  ))
names(optim_param) <- c(paste(rep(names(designTab[,1:length(c(d,m))]),each = 2),c("L","U"),sep=""),"cover")

alpha <- 0.05
nmax_est_param <- length(d)+length(m)
for (i in 1:nrow(designTab)){
  
 index <- which(designTab[i,1:nmax_est_param] %>% as.matrix)
 
 if (length(index) != 0)
 {
 bonf <- length(index)
 
 param <- apply(bootSample, 2, quantile, probs = c(alpha/(2*bonf),1-alpha/(2*bonf))) 
 optim_param[i, c(2*index-1, 2*index)] <- as.numeric(c(param[1, index], param[2, index]))
# optim_param[i, "cover"] = all((param[1,] <= truth & truth <=  param[2,])[,index])
 optim_param[i, "cover"] = all((param[1,] <= truth & truth <=  param[2,])[,index])
  if (!(4 %in% index))  optim_param[i, c("m0L","m0U")] <- c(0,1)
 }else{
   optim_param[i, c("m0L","m0U")] <- c(0,1)
   optim_param[i, "cover"] = TRUE
 }
}

save(bootSample,optim_param,designTab,truth, file = "data/log_Hill_data_optimparam.Rdata")

I = unique(designTab[,1:nmax_est_param]) %>% row.names %>% as.numeric()
write(I, file = "data/log_Hill_data_rows_to_run.txt",sep = "\n")
