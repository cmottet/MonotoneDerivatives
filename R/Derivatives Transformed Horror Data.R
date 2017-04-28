library(plyr)
library(reshape)
library(GLP)
library(magrittr)
library(ggplot2)
library(DistributionPty)
# Load Simulated Data
load("data/Hill_Horror_dist.Rdata")

Tsample <- sort(log(sample))
a <- seq(Tsample[3],Tsample[length(Tsample)-3],length = 100)

###
### Compute the 95% confidence intervals
### for the parameters eta, beta, nu over different 
### values of a
###
# CI <- getCIMomentAndDerivatives(Tsample, 
#                                 a,
#                                 m = 0:4,
#                                 d = 1:3,
#                                 bootSample = TRUE,
#                                 mc.cores = 4)  # Increase the number of cores to use parallel computing

load("data/horrorCI.RData")
save(CI,file = "data/horrorCI.RData")

###
### Transform the CI's in a data frame format
###
library(plyr)
library(dplyr)
dataPlot <- NULL
for (i in 1:length(CI))
{
  bootSample <- CI[[i]]$bootSample
  
  newDataPlot <- data.frame(a =  CI[[i]]$a,
                            parameter = rep(c("First order density derivative", "Second order density derivative", "Density function"),3), 
                            value =  as.numeric(c(CI[[i]]$hyperrectangle[1,1:3], CI[[i]]$hyperrectangle[2,1:3], as.numeric(apply(bootSample[,1:3],2,mean)))), 
                            group = rep(c("lB", "uB", "Fhat"),each  = 3),
                            type = c(rep("Boostrap 95% CI", 6),rep("Boostraped estimated function", 3) ))
  dataPlot <- rbind(dataPlot, newDataPlot)
}


bitmap("pics/EstimationLogHorror.png",res = 300, width = 5,height = 5)
ggplot(dataPlot, aes(x = a, y  = value, group = group)) + 
  geom_line(aes(linetype = type)) + 
  geom_vline(xintercept = qlhorror(0.8), colour = "grey")+
  labs(y = "", linetype = "", x = "") + 
  facet_wrap(~parameter, ncol = 2, scales = "free") + 
  theme(legend.position = c(9.95/10, 2/8), legend.justification = c(1, 0)) 
dev.off()























reshapeBoundsAddTrueValues <- function(data){
  
  output <- NULL
  for (i in 1:length(data)){
    
    df <- data[[i]]$hyperrectangle
    derivatives <- df %>% dplyr::select(grep("d",names(df)))
    moments <- df %>% dplyr::select(grep("m",names(df)))
    
   ciDerivatives <- reshape(derivatives, 
            direction = "long",
            varying = list(1:ncol(derivatives)),
            v.names = "value", 
            ids = c("CI", "CI"),
            idvar = "type",
            timevar = "order", 
            new.row.names = 1:prod(dim(derivatives))) %>%  cbind(data.frame(param = "d",a = data[[i]]$a)) 
  
    trueDerivatives <- ddply( expand.grid(a =a,order =d), 
                              .(a,order), 
                              .fun = function(x) data.frame(value = Dlhorror(x$a,x$order), 
                                                            type = "Truth",
                                                            param = "d"))
    
    
    ciMoments <- reshape(moments, 
                             direction = "long",
                             varying = list(1:ncol(moments)),
                             v.names = "value", 
                             ids = c("CI", "CI"),
                             idvar = "type",
                             timevar = "order",
                             new.row.names = 1:prod(dim(moments))) %>%  cbind(data.frame(param = "m",a = data[[i]]$a)) 
    
    trueMoments <- ddply( expand.grid(a =a,order =m), 
                              .(a,order), 
                              .fun = function(x) data.frame(value = UpperTrunclHorrorMoment(x$a,x$order), 
                                                            type = "Truth",
                                                            param = "m"))
    
    
    
    output <- output %>% rbind(ciDerivatives) %>% rbind(trueDerivatives[names(ciDerivatives)]) %>%
      rbind(ciMoments) %>% rbind(trueMoments[names(ciMoments)])
  }
  row.names(output) = 1:nrow(output)
  return(output)
}


# Add the true values of m and d to the data
reshapedDataAndTruth <- reshapeBoundsAddTrueValues(Data)
ggplot(subset(reshapedDataAndTruth, param == "d"), aes(x = a, y = value)) + 
  geom_line(aes(linetype = type)) +
  facet_grid(order~., scale = "free")


## Derivatives
positive = TRUE
args = list(eval.points = a, positive = positive, func = kcde)
tmp  = CI.bootstrap(Tsample,f1,args,bonferroni = 1,seed = 100)
DataNu = data.frame(a =a, lower = tmp$CI$lower,pred = tmp$Fhat, upper =tmp$CI$upper, Truth = plhorror(a),order = 0)

for (order in 1:3){
  args = list(eval.points = a, positive = FALSE,deriv.order = order-1, func = kdde)
  tmp     = CI.bootstrap(Tsample,f2,args,bonferroni = 1,seed = 100)
  
  ftruth = switch(order,
                  "1" = dlhorror,
                  "2" = ddlhorror,
                  "3" = d2dlhorror)
  DataNutmp = data.frame(a =a,lower = tmp$CI$lower,pred = tmp$Fhat, upper = tmp$CI$upper, Truth = ftruth(a),order = order)
  DataNu = rbind(DataNu,DataNutmp )
}

## Moments
DataMu = NULL
for (order in 1:4){
  fboot =  eval(substitute(function(sample,args) matrix(mean(sample^i*(sample>=args))),list(i=order)))
  tmp = sapply(a,CI.bootstrap,data = Tsample,fboot = fboot,bonferroni = 1,seed = 100)
  DataMutmp  = data.frame(a = a,
                          lower = matrix(unlist(tmp[2,]),ncol=2,byrow = TRUE)[,1],
                          pred = unlist(tmp[1,]), 
                          upper =  matrix(unlist(tmp[2,]),ncol=2,byrow = TRUE)[,2], 
                          Truth = sapply(a,UpperTrunclHorrorMoment,order = order),
                          order = order)
  
  DataMu = rbind(DataMu,DataMutmp )
}

Data = rbind(transform(DataNu,param = "Nu"),
             transform(DataMu,param = "Mu"))
##
## Plot
##, y = value,group = variable)) + 
plot_param_est = function(Data,labelwrap){
  ggplot(data = Data, aes(x = a))+
    geom_ribbon(aes(ymin = lower,ymax = upper),fill = "grey70") + 
    geom_line(aes(y  = Truth))+
    facet_wrap(~order, ncol = 2,nrow = 2, scales = "free_y",
               labeller = as_labeller(labelwrap))+
    geom_vline(xintercept = quantile(Tsample,0.8),linetype = "dashed") + 
    scale_x_continuous(breaks = c(1:3,quantile(Tsample,0.8)), 
                       labels = c(1:3,expression("q"[80])))+
    labs(y = "") +
    theme(axis.text.x = element_text(size = rel(1.2)))
}

plotNu = plot_param_est(DataNu, 
                        c("0" = "Distribution function",
                          "1" = "Density function",
                          "2" = "Density derivative function", 
                          "3" = "Second order density derivative function"))
output.file = file.path(out_pic,"Paper Plots/Derivatives_Transformed_Horror_Data.png")
ggsave(output.file,plot = plotNu, width =10,height = 5) # Save plot in png


plotMu = plot_param_est(DataMu, 
                        c("1" = "1st upper truncated moment",
                          "2" = "2nd upper truncated moment",
                          "3" = "3nd upper truncated moment", 
                          "4" = "4th upper truncated moment"))
output.file = file.path(out_pic,"Paper Plots/Moments_Transformed_Horror_Data.png")
ggsave(output.file,plot = plotMu, width =10,height = 5) # Save plot in png



