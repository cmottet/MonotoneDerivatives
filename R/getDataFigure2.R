remove(list = ls())

# Load Simulated Data
load("data/log_Hill_Horror_dist.Rdata")

sample <- sort(sample)
a <- seq(sample[3],sample[length(sample)-3],length = 100)

###
### Compute the 95% confidence intervals for F(a), F'(a), and F''(a)
### and for different values of a
###
CI <- GLP::getCIMomentAndDerivatives(sample,
                                a,
                                m = 0:4,
                                d = 1:3,
                                bootSample = TRUE,
                                mc.cores = 4)  # Increase the number of cores to use parallel computing

save(CI,file = "data/log_Hill_Horror_CI.RData")
