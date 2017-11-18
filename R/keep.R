# Table containing all the possible scenarios to run the optimization program over
# For example, the first row of designTab indicates that the first optimization
# will be performed with no constraints on the derivatives but it will have a bound
# on the tail distribution and on the first three moments. We will assume no 
# monotonicity for F, i.e. D=0, and the optimization will be performed at the 90
# percentile (P = 0.9).
# Direction = 1 stand for maximization.

# D needs to be larger than the largest derivative,
# we drop the rows that do not meet this constraint
dMax <- apply(designTab[,1:length(d)],1,function(x)max(d[x]))
dMax[dMax == -Inf] <- 0
designTab$dMax <- dMax
designTab <- subset(designTab, dMax <= D) ; rownames(designTab) <- 1:nrow(designTab)

#P <- 1- c(10^-1,10^-2,10^-3,10^-4)
#D <- 0:5
#q <- qlhorror(P)
#direction = 1 # set to c(-1,1) if need to maximize and minimize objective function


