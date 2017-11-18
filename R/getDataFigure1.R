remove(list = ls())

library(DistributionPty)

###
### 0. Set the seed, the number of observations per sample, the number of sample
###
n = 500
seed = 250

###
### Simulate Log Hill Horror Plot distribution <=> F(x) = 1-exp^{-x}/x 
###
set.seed(seed)
sample = rlhorror(n) ; save(sample,file = file.path("data/log_Hill_Horror_dist.Rdata"))

