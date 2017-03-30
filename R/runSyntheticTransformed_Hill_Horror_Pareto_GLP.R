rm(list=ls())

myArgs = commandArgs()
i = as.numeric(myArgs[length(myArgs)]) ; print(i) #n=1

###
### 0. Load files, set seed and define output file
###

# 0.1 Load samples
source("/project/simulate/cmottet/Regularized2/programs/R/init.R") 
input.file  = file.path(input_dir,"data/Hill_Horror_dist.Rdata")
output.file = file.path(out_data,"TransformedHorrorDist_TailProb.Rdata")
load(input.file)
Tsample = log(sample)

if (i == 1)
{
  ###
  ### 1. Set the parameter xf, a and the seed
  ###
  
  a   = as.numeric(quantile(Tsample,0.8)) # Threshold a for GLP
  seed  = 100
  
  ##
  ## 4. True values of all the parameters at a
  ##
  
  nu0 = 1-plhorror(a)
  nu1    = dlhorror(a)
  nu2    = -ddlhorror(a)
  nu3    = d2dlhorror(a)
  mu1 = UpperTrunclHorrorMoment(a,1)
  mu2 = UpperTrunclHorrorMoment(a,2)
  mu3 = UpperTrunclHorrorMoment(a,3)
  mu4 = UpperTrunclHorrorMoment(a,4)
  
  TruthL  = data.frame(nu3,nu2,nu1,nu0)
  TruthM = data.frame(mu1,mu2,mu3,mu4)
  
  ###
  ### 2. Set the parameter to vary and create a design of experiment
  ### Interval over which to maximize
  
  p      = 1- c(10^-1,10^-2,10^-3,10^-4) 
  q      =log(qhorror(p))
  orderM = c(1,2,3,4)
  orderL = 3:0
  
  
  designTab = desnum(fac.design(nlevels=rep(2,length(orderM) + length(orderL) ),randomize = FALSE))
  designTab = designTab==1 
  #   keep = which(rowSums(designTab[,length(orderL):ncol(designTab)]) >= 1)
  #   designTab  = designTab[keep,]
  
  designTabL = designTab[,1:length(orderL)]
  designTabM = designTab[, length(orderL) + (1:length(orderM))]
  
  
  ###
  ### 3. Initialize output
  ###
  names = c(paste("nu",orderL,sep=""), paste("mu",orderM,sep=""), 
            "ParamCover","Primal_uBCover","elapsed.time","CPU.time",
            "primal_uB","dual_lB","status","eps","seed","xf","lastx"
  )
  tmp = data.frame(cbind(designTab,matrix(NA,ncol = length(names)-length(orderM) - length(orderL),nrow = nrow(designTab))))
  names(tmp) = names
  
  optim = rep(list(tmp),length(q))
  
  tmp = list( p = list(NULL), x= list(NULL), r = list(NULL), dualm = list(NULL))
  optim.outDist = rep(list(tmp),length(q))
  
  optim.param = data.frame(matrix(NA,ncol = 2*length(c(orderM,orderL))+1,nrow = nrow(designTab)))
  names(optim.param) = c(paste(rep(c(paste("nu",3:0,sep=""), paste("mu",orderM,sep="")),each = 2),c("L","U"),sep=""),"cover")
  
  const.fun = obj.fun = list(NULL)
  
} else load(output.file)


ind_l = unlist(optim[[1]][i,1:length(orderL)])
ind_m = unlist(optim[[1]][i,length(orderL) + (1:length(orderM))])

l = orderL[ind_l] 
m = orderM[ind_m]

if (length(l)!= 0) L = max(l) else L = 0


##
## Create the set of function constraints and objective function
##

const.fun[[i]] = rep(buildMomentDerConstFunc(L,l,m),each = 2)
obj.fun[[i]]   = eval(substitute (function(x,param.obj.fun)   # H for the optim with convexity constraint
{
  a = param.obj.fun$a
  c = param.obj.fun$c
  
  output = 1/factorial(L)*(x - max(c-a,0))^L*( max(c-a,0) <= x) 
  return(output)
}  ,
list(L=L)) ) 

###
### The values on the RHS of the constraints
###
tmp = getOptimParam(sample = Tsample,
                    a = a,
                    orderM = m,
                    orderL = l,
                    TruthM = TruthM[ind_m],
                    TruthL = TruthL[ind_l],
                    seed = seed )

# tmp = optim.param[i,which(!is.na(optim.param[i,1:15]))]

optim.param[i,which(names(optim.param)%in% names(tmp))] = tmp

const.rhs = tmp
if ("cover" %in% names(const.rhs)) const.rhs = subset(const.rhs,select = -cover)


###
## Sense of the Inequality 
##
const.dir = rep(c(">=", "<="), length(const.rhs)/2) 


##
## Create the set of lambdas
##
lambda.const = rep(0,length(const.rhs)) 
lambda.const[length(const.rhs) - c(1,0)] = 1

if (length(m) > 0){
  lambda.obj = 0
} else lambda.obj = 1


##
## Compute xf
##
# xf = 10E10
# if (length(m)!= 0 ) d = sort(m)[length(m)] + L  else d = L
# if (d != 0) xf = xf^(1/d)
xf = log(qhorror(1-1e-16))
##
## Compute the initial feasible solution
##

initBFS  =  GLP.phase1(const.fun[[i]],
                       const.rhs,
                       const.dir,
                       param.const.fun = a,
                       xf =  xf, IterMax = 200)# xf = min(xf,400)


for (j in 1:length(q))
{
     
  param.obj.fun = list(a = a, c= q[j])
  obj.func.indic = sum(unlist(optim[[j]][i,1:3])) == 0
  
  # Compute optimal upper bound and measure the time needed to do so
  start = proc.time()    
  out = GLP.phase2(initBFS,
                   obj.fun[[i]], 
                   const.fun[[i]],
                   const.rhs,
                   const.dir,
                   lambda.const,
                   lambda.obj,
                   param.const.fun = a,
                   param.obj.fun,
                   xf = xf,
                   IterMax = 200,
                   obj.func.indic = obj.func.indic) 
  end = proc.time()
  
  # Record results
  optim[[j]]$elapsed.time[i] = (end - start)["elapsed"]
  optim[[j]]$CPU.time[i]     = (end - start)["user.self"]
  
  optim[[j]]$primal_uB[i] = out$primal_uB
  optim[[j]]$dual_lB[i]   = out$dual_lB
  optim[[j]]$status[i]    = out$status
  optim[[j]]$eps[i]       = out$eps         
  optim[[j]]$seed[i]      = seed   
  optim[[j]]$ParamCover[i] = optim.param$cover[i]
  optim[[j]]$xf[i]        = xf
  optim[[j]]$Primal_uBCover[i] = out$primal_uB >=1-p[j]
  optim[[j]]$lastx[i] = out$lastx
  
  optim.outDist[[j]]$x[[i]] = out$x
  optim.outDist[[j]]$p[[i]] = out$p
  optim.outDist[[j]]$r[[i]] = out$r
  optim.outDist[[j]]$dualm[[i]] = out$lpdual
  
  
  ###
  ### Save the results and display    
  ###
  print(c(l,m,out$dual_lB,out$primal_uB,1-p[j]))
  save(optim,
       optim.param, 
       optim.outDist,
       q,p,a,
       TruthM,
       TruthL,
       orderM,
       orderL,
       seed,
       obj.fun,
       const.fun,
       file= output.file)
}

# lpdual = optim.outDist[[j]]$dualm[[i]]
# 
# # inOpt= GenSA(par       = 0,
# #       fn        = finOpt,# GenSA for GLOBAL max
# #       lower     = 0,
# #       upper     = optim[[j]]$xf[i],
# #       phase     = 2,
# #       lpdual    = lpdual,
# #       const.fun = const.fun[[i]],
# #       param.const.fun = a,                  
# #       obj.fun   = obj.fun[[i]],
# #       param.obj.fun   = param.obj.fun)
# 
# # inOpt$value
# xf = optim[[j]]$xf[[i]]
# lastx = optim[[j]]$lastx[i]
# eps = optim[[j]]$eps[i]
# 
# library(emdbook)
# 
# xmin = 906.79
# xmax = 906.8
# x= lseq(xmin,xmax,1e4)
# y = sapply(x,finOpt,
#            lpdual = lpdual,
#            phase = 2,
#            const.fun = const.fun[[i]],
#            param.const.fun = a,
#            obj.fun = obj.fun[[i]],
#            param.obj = param.obj.fun)
# 
# Data = data.frame(x=x,y=y)
# plot = ggplot(Data,aes(x=x,y=y)) + geom_line()
# plot