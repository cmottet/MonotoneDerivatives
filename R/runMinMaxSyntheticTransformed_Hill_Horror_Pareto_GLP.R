rm(list=ls())


myArgs = commandArgs()
i = as.numeric(myArgs[length(myArgs)]) ; print(i) #n=1

###
### 0. Load files, set seed and define output file
###
library(GLP)
library(DistributionPty)
library(dplyr)

setwd("/home/ec2-user/MonotoneDerivatives")
source("R/extras.R")
load("data/log_Hill_data_optimparam.Rdata") # Load optim_param, design_tab
load("data/Hill_Horror_dist.Rdata")

Tsample = log(sample)
# load(output.file)
# tmp =which(apply(optim[,13:17],1,function(x)sum(!is.na(x))>0 ))
# length(which(apply(optim[,13:17],1,function(x)sum(!is.na(x))>0 )))

a <- as.numeric(quantile(Tsample,0.8)) # Threshold a for GLP
m <- 0:4
d <- 3:1
D <- 0:5
P <- 1- c(10^-1,10^-2,10^-3,10^-4)

nmax_est_param <- length(m) + length(d)

if (i == 1)
{
  ### 1. Initialize output
  names = c(names(designTab),
            "ParamCover","boundCover","elapsed.time","CPU.time",
            "uB","lB","status","eps","seed","xf","lastx")
  
  constFun = objFun = list()
  optim_outDist = list( p = list(), x= list(), r = list(), dualm = list())
  optim <- data.frame(matrix(NA,ncol = length(names),nrow = nrow(designTab))) 
  optim[,1:ncol(designTab)] = designTab
  names(optim) <- names
  optim$xf <- log(qhorror(1-1e-16))
  
} else  load("data/MinMaxTransformedHorrorDist_TailProb.Rdata")


# Sense of the Inequality and RHS of constraints
constRHS = optim_param[i,!is.na(optim_param[i,])] %>% select(-cover) 
constDir = rep(c(">=", "<="), length(constRHS)/2)

# Find all the rows in the DoE that have the parameters associated with i
param_est <- designTab[i,1:nmax_est_param]
K <- which(apply(designTab[,1:nmax_est_param],1, function(x) all(x== param_est)))

initBFS <- rep(list(NULL),length(D))
for (k in K)
{
  ind_d <- unlist(param_est[1:length(d)])
  ind_m <- unlist(param_est[(length(d)+1):nmax_est_param])
  
  # Create the set of function constraints and objective function (change with D)
  constFun[[k]] <- rep(buildMomentDerConstFunc(designTab$D[k],d[ind_d],m[ind_m]),each = 2)
  
  # Create Lambda const
  constLambda = rep(0,length(constRHS))
  constLambda[length(constRHS) - c(1,0)] = 1
  
  # Compute the initial feasible solution
  if (is.null(initBFS[[optim$D[k]+1]]))
    initBFS[[optim$D[k]+1]] =  GLPPhase1(constFun[[k]],
                                         constRHS,
                                         constDir,
                                         paramConsFun = a,
                                         xf =  optim$xf[k],
                                         IterMax = 200)# xf = min(xf,400)
  
  # Objective function
  objFun[[k]]   = eval(substitute (function(x,paramObjFun)   # H for the optim with convexity constraint
  {
    a = paramObjFun$a
    c = paramObjFun$c
    
    output = direction/factorial(D)*(x - max(c-a,0))^D*( max(c-a,0) <= x)
    return(output)
  }  ,
  list(D= optim$D[k], direction = optim$direction[k])
  ))
  
  # Lambda of the objective function
  objLambda <- if (length(m[ind_m]) > 0)  0 else optim$direction[k]
  
  paramObjFun <- list(a = a, c= qlhorror(optim$P[k]))
  objFuncIndic <- optim$D[k] == 0
  
  gamma <- NULL
  if (optim$D[k] %in% d[ind_d]) gamma <- optim_param[i,paste0("d",optim$D[k],"U")]
  if (optim$D[k] == 0) gamma <- optim_param[i,paste0("m0U")]
  
  # Compute optimal upper bound and measure the time needed to do so
  start <- proc.time()
  out <- GLPPhase2(initBFS = initBFS[[optim$D[k]+1]],
                   objFun = objFun[[k]],
                   constFun = constFun[[k]],
                   constRHS,
                   constDir,
                   constLambda,
                   objLambda,
                   paramConsFun = a,
                   paramObjFun,
                   gamma = gamma,
                   xf = optim$xf[k],
                   IterMax = 200,
                   objFuncIndic = objFuncIndic,
                   factor = optim$direction[k])
  end <- proc.time()
  
  # Record results
  optim$elapsed.time[k] <- (end - start)["elapsed"]
  optim$CPU.time[k]     <- (end - start)["user.self"]
  optim$uB[k]         <- out$uB
  optim$lB[k]         <- out$lB
  optim$status[k]     <- out$status
  optim$eps[k]        <- out$eps
  optim$ParamCover[k] <- optim_param$cover[i]
  optim$lastx[k]      <- out$lastx
  
  optim$boundCover[k] <- if (optim$direction[k] == 1)  out$lB >= 1-optim$P[k] else out$uB <= 1-optim$P[k]
  
  optim_outDist$x[[k]] <- out$x
  optim_outDist$p[[k]] <- out$p
  optim_outDist$r[[k]] <- out$r
  optim_outDist$dualm[[k]] <- out$lpdual
  
  # Save the results and display
  sum_print <- matrix(c(optim$direction[k],optim$D[k], d[ind_d],m[ind_m],signif(out$lB,5),signif(out$uB,5),1-optim$P[k] ),nrow = 1)
  colnames(sum_print) <- c("Dir","D",names(ind_d[ind_d]),names(ind_m[ind_m]),"lB","uB", "P" )
  print(sum_print)
  
  save(optim,optim_param,optim_outDist,
       q,P,a,D,
       m,d,
       objFun,constFun,
       file= "data/MinMaxTransformedHorrorDist_TailProb.Rdata")
}
