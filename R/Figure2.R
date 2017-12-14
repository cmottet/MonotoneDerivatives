remove(list = ls())
library(ggplot2)

# Load CI's
load("data/log_Hill_Horror_CI.RData")

#
# Transform the CI's in a data frame format
#
dataPlot <- NULL
for (i in 1:length(CI))
{
  bootSample <- CI[[i]]$bootSample
  
  newDataPlot <- data.frame(a =  CI[[i]]$a,
                            parameter = rep(c("Second order density derivative", "First order density derivative", "Density function"),3), 
                            value =  as.numeric(c(CI[[i]]$hyperrectangle[1,1:3], CI[[i]]$hyperrectangle[2,1:3], as.numeric(apply(bootSample[,1:3],2,mean)))), 
                            group = rep(c("lB", "uB", "Fhat"),each  = 3),
                            type = c(rep("Bootstrap 95% CI", 6),rep("Bootstrapped estimated function", 3) ))
  dataPlot <- rbind(dataPlot, newDataPlot)
}

#
# Plot the CI's
#
plot <- ggplot(dataPlot, aes(x = a, y  = value, group = group)) + 
  geom_line(aes(linetype = type)) + 
  geom_vline(xintercept = DistributionPty::qlhorror(0.8), colour = "grey")+
  labs(y = "", linetype = "", x = "") + 
  facet_wrap(~parameter, ncol = 2, scales = "free") + 
  theme(legend.position = c(9.95/10, 2/8), legend.justification = c(1, 0)) 

plot 

ggsave(plot, file = "pics/Figure2.pdf", width = 5,height = 5,dpi=300)



