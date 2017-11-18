rm(list=ls())
library(dplyr)
source("R/extras.R")


myArgs <- commandArgs()
i <- as.numeric(myArgs[length(myArgs)]) ; print(i) #n=1

###
### 0. Load files, set seed and define output file
###

load("data/log_Hill_data_optimparam.Rdata") # Load optim_param
load("data/log_Hill_Horror_dist.Rdata")

a <- as.numeric(quantile(sample,0.8)) # Threshold a for GLP
m <- 0:3
d <- 3:1
P <- 1- c(10^-1,10^-2,10^-3,10^-4)

nmax_est_param <- length(m) + length(d)

#
# 1. Initialize output
#
if (i == 1)
{
  names = c(names(designTab), "D", "P", "direction",
            "ParamCover","boundCover","elapsed.time","CPU.time",
            "uB","lB","status","eps","C","lastx")
  
  optim <- data.frame(matrix(NA,ncol = length(names), nrow = 1)) 
  names(optim) <- names
  
  constFun = objFun = optim_outDist = list()
  
} else  load("data/MaxLogdHorrorDist_TailProb.Rdata")

# Find all the rows in the DoE that have the parameters associated with i
param_est <- designTab[i,]

# Index of estimated derivatives and moments
ind_d <- unlist(param_est[1:length(d)])
ind_m <- unlist(param_est[(length(d)+1):nmax_est_param])

# Obtain the set J1 and J2
J1 <- union(0, m[ind_m]) 
J2 <- d[ind_d]

# Set the value of J
J <- max(J2)
J <- if (J == -Inf)  0 else J

C <- DistributionPty::qlhorror(1-1e-16) # Upper bound of the compact support of the distributions to investigate from
direction <- 1 # Maximize optimization program

# Get bounds on constraints inequalities
point_estimates <- optim_param[i,!is.na(optim_param[i,])] %>% select(-cover) 
constRHS <- NULL

for (j in J2){
  newRHS <- point_estimates[paste0("d",j,c("L","U"))] %>% as.numeric %>% `*`((-1)^(j+1)) %>% sort
  constRHS <- c(constRHS, newRHS)
}  

for (j in J1){
  newRHS <- point_estimates[paste0("m",j,c("L","U"))] %>% as.numeric
  constRHS <- c(constRHS, newRHS)
}

# Get the direction of constraints ineqaulities
constDir = rep(c(">=", "<="), length(constRHS)/2)

# Create Lambda const
constLambda = rep(0,length(constRHS))
constLambda[length(constRHS) - c(1,0)] = 1

# Lambda of the objective function
objLambda <- if (setequal(J1,0))  1 else 0

# Get the upper bound on the distribution function
gamma <- constRHS[2] 

k <- nrow(optim)
for(D in J:5)
{
  # Is the objective function an indicator?
  objFuncIndic <- D == 0
  
  # Create the set of function constraints and objective function (change with D)
  new_constFun <- rep(buildMomentDerConstFunc(D,J1,J2), each = 2)
  
  # Compute the initial feasible solution
  initBFS <- GLP:: GLPPhase1(new_constFun,
                             constRHS,
                             constDir,
                             paramConsFun = a,
                             xf =  C,
                             IterMax = 200)
  
  for (p in P)
  {
    # Parameters of the objective function
    paramObjFun <- list(a = a, c = DistributionPty::qlhorror(p))
    
    # Objective function
    objFun[[k]]  <- function(x, paramObjFun) 
    {
      a <- paramObjFun$a
      c <- paramObjFun$c
      
      output <-  if (max(c-a,0) < x) direction*x^(J-D)/factorial(D)*(x - max(c-a,0))^D else 0
      return(output)
    } 
    
    constFun[[k]] <- new_constFun
    
    # Compute optimal upper bound and measure the time needed to do so
    start <- proc.time()
    out <- GLP::GLPPhase2(initBFS = initBFS,
                          objFun = objFun[[k]] ,
                          constFun = constFun[[k]] ,
                          constRHS,
                          constDir,
                          constLambda,
                          objLambda,
                          paramConsFun = a,
                          paramObjFun,
                          gamma = gamma,
                          xf = C,
                          IterMax = 200,
                          objFuncIndic = objFuncIndic,
                          factor = direction)
    end <- proc.time()
    
    # Record results
    if (k!=1) optim <- rbind(optim, NA) # Add a row if k is not 1
    
    optim[k,1:nmax_est_param] <- designTab[i, ]
    optim$D[k] <- D
    optim$P[k] <- p
    optim$direction[k] <- direction
    optim$ParamCover[k] <- optim_param$cover[i]
    optim$elapsed.time[k] <- (end - start)["elapsed"]
    optim$CPU.time[k]     <- (end - start)["user.self"]
    optim$uB[k]         <- out$uB
    optim$lB[k]         <- out$lB
    optim$status[k]     <- out$status
    optim$eps[k]        <- out$eps
    optim$lastx[k]      <- out$lastx
    optim$C[k]  <- C
    optim$boundCover[k] <- if (optim$direction[k] == 1)  out$lB >= 1-optim$P[k] 
    
    optim_outDist[[k]] <- list(x = out$x, p = out$p,  s = out$r, lpdual = out$lpdual)
    
    # Display and save the results 
    sum_print <- matrix(c(direction,D,d[ind_d],m[ind_m],signif(out$lB,5),signif(out$uB,5),1-p),nrow = 1)
    colnames(sum_print) <- c("Dir","D",names(ind_d[ind_d]),names(ind_m[ind_m]),"lB","uB", "P" )
    print(sum_print)
    print(" ")
    
    save(optim,optim_param,optim_outDist,
         P,a,m,d,truth,
         objFun,constFun,
         file = "data/MaxLogdHorrorDist_TailProb.Rdata")
    
    k <- k+1
  }
}
