# !!!!!! Run this script from the command line by typing !!!!!
# ./shell/RjobRun.sh
# from the home directory of this repositary

rm(list=ls())
library(dplyr)

myArgs <- commandArgs()
i <- as.numeric(myArgs[length(myArgs)]) ; print(i)

###
### 0. Load files, set seed and define output file
###

load("data/log_Hill_data_optimparam.Rdata") # Obtain by running getOptimParam.R
load("data/log_Hill_Horror_dist.Rdata") # Obtain by running getDataFigure1.R

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
  names = c(names(designTab), "D", "P",
            "ParamCover","boundCover","elapsed.time","CPU.time", "bound","status","eps","C","lastx")

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

k <- nrow(optim)
for(D in J:5)
{
  # Create the set of function constraints and objective function (change with D)
  new_constFun <- rep(GLP::buildMomentDerivativeConstFunc(D,J1,J2), each = 2)

  # Compute the initial feasible solution
  initBFS <- GLP:: phase1(new_constFun,
                          constRHS,
                          constDir,
                          constLambda,
                          C =  C,
                          IterMax = 200)

  for (p in P)
  {
    # Objective function
    objFun[[k]]  <- function(x, paramObjFun)
    {
      c <- DistributionPty::qlhorror(p)

      output <-  if (max(c-a,0) < x) x^(J-D)/factorial(D)*(x - max(c-a,0))^D else 0
      return(output)
    }

    constFun[[k]] <- new_constFun

    # Compute optimal upper bound and measure the time needed to do so
    start <- proc.time()
    out <- GLP::phase2(initBFS = initBFS,
                       objFun = objFun[[k]] ,
                       constFun = constFun[[k]] ,
                       constRHS,
                       constDir,
                       constLambda,
                       objLambda,
                       C = C,
                       IterMax = 200,
                       err = 1e-7)
    end <- proc.time()

    # Record results
    if (k!=1) optim <- rbind(optim, NA) # Add a row if k is not 1

    optim[k,1:nmax_est_param] <- designTab[i, ]
    optim$D[k] <- D
    optim$P[k] <- p
    optim$ParamCover[k] <- optim_param$cover[i]
    optim$elapsed.time[k] <- (end - start)["elapsed"]
    optim$CPU.time[k]     <- (end - start)["user.self"]
    optim$bound[k]      <- out$bound
    optim$status[k]     <- out$status
    optim$eps[k]        <- out$eps
    optim$lastx[k]      <- out$lastx
    optim$C[k]  <- C
    optim$boundCover[k] <- out$bound >= 1-optim$P[k]

    optim_outDist[[k]] <- list(x = out$x, p = out$p,  s = out$s, lpdual = out$lpdual)

    # Display and save the results
    sum_print <- matrix(c(D,d[ind_d],m[ind_m],signif(out$bound,5),1-p),nrow = 1)
    colnames(sum_print) <- c("D",names(ind_d[ind_d]),names(ind_m[ind_m]),"bound", "P" )
    print(sum_print)
    print(" ")

    save(optim,optim_param,optim_outDist,
         P,a,m,d,truth,
         objFun,constFun,
         file = "data/MaxLogdHorrorDist_TailProb.Rdata")

    k <- k+1
  }
}

