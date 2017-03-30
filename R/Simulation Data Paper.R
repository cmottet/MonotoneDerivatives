source("/project/simulate/cmottet/Regularized2/programs/R/init.R") 

###
### 0. Set the seed, the number of observations per sample, the number of sample
###
nobs = 200
nsam = 100
seed = 100

###
### Simulate a Pareto sample with shape = 2, location = 1
###

dist = "Pareto"
alpha = 2 ; xi = 1/alpha ; beta = 1/alpha; loc = 1 # Location of the Pareto
set.seed(seed) ; samples = matrix(loc + rGPD(nsam*nobs,xi = xi,beta = beta),nrow = nsam)
save(samples, xi,beta,loc,dist,file = file.path(input_dir,"data/Pareto.Rdata"))

# pGPD(max(apply(samples,1,quantile,prob = 0.8)) - loc,xi = xi,beta=beta) # 0.8577

###
### Simulate an Exponential Sample with rate = 2
###

dist = "Exponential"
rate = 2 ; 
set.seed(seed) ; samples = matrix(rexp(nsam*nobs,rate),nrow = nsam)
save(samples,rate,dist,file = file.path(input_dir,"data/Exponential.Rdata"))
# pexp(max(apply(samples,1,quantile,prob = 0.8)),rate) #0.87

###
### Simulate a Log-Normal Sample with mean 0 and meanlog 0.5
###

dist = "Log-Normal"
meanlog = 0 ; sdlog = 0.5 
set.seed(seed) ; samples = matrix(rlnorm(nsam*nobs,meanlog,sdlog),nrow = nsam)
save(samples,meanlog,sdlog,dist,file = file.path(input_dir,"data/LogNormal.Rdata"))
# plnorm(max(apply(samples,1,quantile,prob = 0.8)),meanlog,sdlog) #0.857


###
### Simulate Hill Horror Plot distribution
###

n = 500
sample = rhorror(n) ; save(sample,file = file.path(input_dir,"data/Hill_Horror_dist.Rdata"))